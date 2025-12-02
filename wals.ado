* Weighted Average Least Squares (WALS) with iid errors
*! Version 3.0
* Date: 2024/10/05
* Authors: De Luca Giuseppe, Jan Magnus, 
* Version 1.0: De Luca G. and Magnus J.R. (2011)  
* Version 2.0: De Luca G. and Magnus J.R. (2016)  
* Version 3.0: De Luca G. and Magnus J.R. (2024)  
*-----------------------------------------------------------------------------------------------------------------------	
* wals
*-----------------------------------------------------------------------------------------------------------------------	
program define wals, eclass sort
	version 14.0, missing
	if replay() {
		if "`e(cmd)'" != "wals" {
			error 301
		}
		syntax [, SAVing(string asis) noHEADer noTABle noFOCus noAUXiliary LEVel(cilevel) cformat(string)]
		local bcsimdata `"`e(bcsimdata)'"'
		cap confirm file `"`bcsimdata'"'
		if _rc {
			di as err "MC replicates of bias-corrected wals estimates not found. You need to refit the model and resave the results" 
			exit _rc
		}
		if `"`saving'"' != "" {
			_savingopt_parse CIMC_save CIMC_savereplace : saving ".dta" `"`saving'"'
			if `"`CIMC_save'"'!="" {
			    if `"`CIMC_savereplace'"'=="" {
				    confirm new file `"`CIMC_save'"'
				}
				if `"`bcsimdata'"' != `"`CIMC_save'"' {
					qui copy `"`bcsimdata'"' `"`CIMC_save'"', `CIMC_savereplace'
					mata: st_global("e(bcsimdata)",`"`CIMC_save'"')
					if `"`bcsimdata'"'==`"`c(tmpdir)'/STWALS_000001.tmp"' cap erase `"`c(tmpdir)'/STWALS_000001.tmp"'
				}
			}	
			if `"`CIMC_save'"'=="" {
				local CIMC_save `"`c(tmpdir)'/STWALS_000001.tmp"'
				local CIMC_savereplace "replace"
			}
		}
		if "`cformat'"!="" local cformat "cformat(`cformat')"
		Di_wals, `header' `table' `focus' `auxiliary' lev(`level') `cformat'
		exit
	}
	if !replay() {
		syntax [anything] [fw aw iw /] [if] [in], [*]
		local cmdline : copy local 0
	}
	Estimate `0'
	ereturn local cmdline `"wals `0'"'
end
*-----------------------------------------------------------------------------------------------------------------------	



*-----------------------------------------------------------------------------------------------------------------------	
* Estimate
*-----------------------------------------------------------------------------------------------------------------------	
program define Estimate, eclass
qui {



*------------------------------------------------------------------------------------------------
* Outcome variable 
*------------------------------------------------------------------------------------------------
gettoken depvar 0 : 0
_fv_check_depvar `depvar'
*------------------------------------------------------------------------------------------------

	

*------------------------------------------------------------------------------------------------
* Syntax 
*------------------------------------------------------------------------------------------------
#delimit;	
syntax [anything] [fw aw iw /] [if] [in] , 											
	[														
	AUXCONstant												
	noCONstant												
	nocollcheck		/* undocumented */
	sigma(real 0)  											
	PRIor(string) 	
	
	QUADMethod(string)					
	QUADNPts(integer 500)		
	QUADExt(string)					
	QUADATol(real 1e-9) 
	QUADRTol(real 1e-7)		
	
	PLUGin(string)											
	PATHTab(string asis)
	
	LEVel(cilevel)
	REPs(numlist max=1 >0 integer)						
	rseed(numlist max=1 >0 integer)						
	SAVing(string asis)										
															
	noHEADer
	noFOCus
	noAUXiliary
	noTABle					
	cformat(string)
	fast
	
	];
#delimit  cr	
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Parse lists of focus and auxiliary regressors
*------------------------------------------------------------------------------------------------
noi wals_parse_regressors `anything'
local foclist0 	`s(foc)'
local auxlist0 	`s(aux)'
if "`:list dups foclist0'"!="" {  
	di as err "focus regressors cannot contain duplicate variables"
	error 198
}
if "`:list dups auxlist0'"!="" {  
	di as err "auxiliary regressors cannot contain duplicate variables"
	error 198
}
if "`:list foclist0 & auxlist0'"!="" {  
	di as err "focus and auxiliary regressors cannot contain duplicate variables"
	error 198
}
if `:list depvar in foclist0' {
	di as err "dependent variable {bf:`depvar'} not allowed as focus regressor"
	error 198
}
if `:list depvar in auxlist0' {
	di as err "dependent variable {bf:`depvar'} not allowed as auxiliary regressor"
	error 198
}
if subinstr("`foclist0' `auxlist0'",".`depvar' ","",.) != "`foclist0' `auxlist0'" {
	di as err "time-series operators of dependent variable {bf:`depvar'} may not be " /*
	*/ 	"included as independent variables"
	error 198
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Mark estimation sample 
*------------------------------------------------------------------------------------------------
marksample touse
markout `touse' `depvar' `foclist0' `auxlist0' `exp'  
count if `touse' 
if r(N)==0 	error 2000 
else 		local N0=r(N)	
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Weights 
*------------------------------------------------------------------------------------------------
if `"`weight'"' != "" & `N0'==1 {
	noi di as error "weights not allowed with one observation"
	exit 198
}
tempvar wgt_var
if `"`weight'"' != "" {
	gen double `wgt_var'=`exp' 					if `touse'
	if `"`weight'"'=="pweight" {
		noi di as err "`weight' not allowed"
		error 198
	}
	else if `"`weight'"'=="aweight" {
		sum `wgt_var' , meanonly
		replace `wgt_var' = `wgt_var'/r(mean)	if `touse'
		local wgt_exp `"[aw=`wgt_var']"'
		local wgt_exp_post `"[aw=`exp']"'		
	}
	else if `"`weight'"'=="fweight" {
		local wgt_exp `"[fw=`wgt_var']"'
		local wgt_exp_post `"[fw=`exp']"'
	}
	else {
		local wgt_exp `"[iw=`wgt_var']"'
		local wgt_exp_post `"[iw=`exp']"'
	}
}
else {
	gen double `wgt_var' = `touse' 				if `touse'
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* wals_options:  prior (default pareto)
*------------------------------------------------------------------------------------------------
local P=lower(trim(`"`prior'"'))
local l = length(`"`P'"')
if `"`P'"'=="" 										local PRI "pareto" 
else if `"`P'"'==substr("laplace"  ,1,max(`l',3))	local PRI "laplace"
else if `"`P'"'==substr("weibull"  ,1,max(`l',3))	local PRI "weibull"
else if `"`P'"'==substr("subbotin" ,1,max(`l',3))	local PRI "subbotin"
else if `"`P'"'==substr("pareto"   ,1,max(`l',3))	local PRI "pareto"
else if `"`P'"'==substr("cauchy"   ,1,max(`l',3))	local PRI "cauchy"
else if `"`P'"'==substr("log"      ,1,max(`l',3))	local PRI "log"
else if `"`P'"'==substr("horseshoe",1,max(`l',3))	local PRI "horseshoe"
else {
	noi di in red "{bf: prior} must be one of the following: laplace, subbotin, weibull, pareto, cauchy, log, horseshoe"
	noi di in red "default prior: pareto"
	error 198
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* wals_options: quadmethod (default gauss - i.e. Gauss-Legendre or Gauss-Laguerre)
*------------------------------------------------------------------------------------------------
local quadm = lower(trim("`quadmethod'"))
if "`quadm'"==""		local quadm="gauss"
if 	"`quadm'"!="gauss" 		& 	///
	"`quadm'"!="adaptive" 	{
	noi di in red "{bf:quadmethod} must be one of the following: gauss, adaptive"
	noi di in red "default quadmethod: gauss"
	error 198
}
if "`quadm'"=="gauss" 	{
	if "`quadext'"!="" {
		cap mata: `quadext'
		if _rc!=0 {
			noi di in red "mata matrix {bf:quadext(`quadext'}) not found"
			error 3499
		}
		else {
			if "`quadnpts'"!="500"	{
				noi di as text "note: {bf:quadnpts} is ineffective with {bf:quadext}"
			}
			tempname wals_GL_PW
			mata: `wals_GL_PW'	=valofexternal(st_local("quadext"))
			mata: st_matrix("r(GL_QP)"		, rows(`wals_GL_PW')	)
			local quadnpts=r(GL_QP)[1,1]
			if "`PRI'"=="horseshoe" {
				mata:  st_local("check_ext_p", strofreal(colmaxabs(`wals_GL_PW'[.,1])<=1))
				if "`check_ext_p'"!="1" {
					noi di as err "Mata matrix {bf:quadext(`quadext')} invalid." /*
						*/ _n "The `PRI' prior requires that abs(max `quadext'[.,1])<=1"
					error 9
				}
				mata:  st_local("check_ext_w", strofreal((colmin(`wals_GL_PW'[.,2])>=0) & (colmax(`wals_GL_PW'[.,2])<1)) )
				if "`check_ext_w'"!="1" {
					noi di as err "Mata matrix {bf:quadext(`quadext')} invalid." /*
						*/ _n "The quadrature weights in `quadext'[.,2] must be in [0,1)"
					error 9
				}
			}
			else {
				mata:  st_local("check_ext_p", strofreal(colmin(`wals_GL_PW'[.,1])>=0))
				if "`check_ext_p'"!="1" {
					noi di as err "Mata matrix {bf:quadext(`quadext')} invalid." /*
						*/ _n "The `PRI' prior requires that (min `quadext'[.,1])>=0"
					error 9
				}
				mata:  st_local("check_ext_w", strofreal((colmin(`wals_GL_PW'[.,2])>=0) & (colmax(`wals_GL_PW'[.,2])<1)) )
				if "`check_ext_w'"!="1" {
					noi di as err "Mata matrix {bf:quadext(`quadext')} invalid." /*
						*/ _n "The quadrature weights in `quadext'[.,2] must be in [0,1)"
					error 9
				}
			}
		}
	}
	if `quadatol'!=1e-9   	noi di as text "note: {bf:quadatol} is ineffective with {bf:quadmethod(gauss)}"
	if `quadrtol'!=1e-7		noi di as text "note: {bf:quadrtol} is ineffective with {bf:quadmethod(gauss)}"
	local COD_quadm=1
}
if "`quadm'"=="adaptive" {		
	if "`quadext'"!=""		noi di as text "note: {bf:quadext} is ineffective with {bf:quadmethod(adaptive)}"
	if "`quadnpts'"!="500"	noi di as text "note: {bf:quadnpts} is ineffective with {bf:quadmethod(adaptive)}"
    local COD_quadm=2
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Setup of prior parameters (based on neutrality and minmax regret values) 
*------------------------------------------------------------------------------------------------
if "`PRI'"=="laplace" {
	tempname a b c prior_par 
	local COD_prior=1
	scalar `c'=log(2)
	if "`quadm'"!="gauss"			noi di as text "note: {bf:quadmethod} is ineffective with {bf:prior(`PRI')}"
	if "`quadext'"!="" 				noi di as text "note: {bf:quadext} is ineffective with {bf:prior(`PRI')}"
	else {
	    if "`quadnpts'"!="500" 		noi di as text "note: {bf:quadnpts} is ineffective with {bf:prior(`PRI')}"
	}
	if `quadatol'!=1e-9 			noi di as text "note: {bf:quadatol} is ineffective with {bf:prior(`PRI')}"
	if `quadrtol'!=1e-7 			noi di as text "note: {bf:quadrtol} is ineffective with {bf:prior(`PRI')}"
	local quadm 	"gauss"
	local quadext 		""
	local quadnpts	=2
	local quadatol	=1e-9
	local quadrtol	=1e-7
	matrix `prior_par'=(`c')
	matrix rown `prior_par'= c
}
if "`PRI'"=="subbotin" {
	tempname a b c prior_par 
	local COD_prior=2
	scalar `b'=0.799442158785679
	scalar `c'=0.93777926867646
	matrix `prior_par'=(`b' \ `c')
	matrix rown `prior_par'= b c
}
if "`PRI'"=="weibull" {
	tempname a b c prior_par 
	local COD_prior=3
	scalar `b'=0.887395795667800
	scalar `c'=log(2)
	matrix `prior_par'=(`b' \ `c')
	matrix rown `prior_par'= b c
}
if "`PRI'"=="pareto" {
	tempname a b c prior_par 
	local COD_prior=4
	scalar `a'=0.086210823981228
	scalar `c'=0.067580105532787
	matrix `prior_par'=(`a' \ `c')
	matrix rown `prior_par'= a c
}
if "`PRI'"=="cauchy" {
	tempname prior_par 
	local COD_prior=5
	matrix `prior_par'=(J(1,1,.))
	matrix rown `prior_par'= none
}
if "`PRI'"=="log" {
	tempname c prior_par 
	local COD_prior=6
	scalar `c'=2.445828316431493
	matrix `prior_par'=(`c')
	matrix rown `prior_par'= c
}
if "`PRI'"=="horseshoe" {
	tempname c prior_par 
	local COD_prior=7
	scalar `c'=1.732872007622451
	matrix `prior_par'=(`c')
	matrix rown `prior_par'= c
}
matrix coln `prior_par'= prior_par
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Load tabulated sampling moments of the posterior mean in the normal location model
*------------------------------------------------------------------------------------------------
if "`fast'"=="" {
	if `COD_prior'==1 	local TABMOM "wals_pm_tsm_laplace.dta"
	if `COD_prior'==2 	local TABMOM "wals_pm_tsm_subbotin.dta"
	if `COD_prior'==3 	local TABMOM "wals_pm_tsm_weibull.dta"
	if `COD_prior'==4 	local TABMOM "wals_pm_tsm_pareto.dta"
	if `COD_prior'==5 	local TABMOM "wals_pm_tsm_cauchy.dta"
	if `COD_prior'==6 	local TABMOM "wals_pm_tsm_log.dta"
	if `COD_prior'==7 	local TABMOM "wals_pm_tsm_horseshoe.dta"
	preserve
	findfile `TABMOM', path(`"`pathtab'"')
	use `"`r(fn)'"', clear
	noi mata: tabmom=st_data(.,"theta PM_BIAS PM_VAR")
	restore
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Options for plug-in estimators (ds and ml) of the sampling moments (bias and variance)
*------------------------------------------------------------------------------------------------
if "`fast'"=="" {
	local BIAS_type = lower(trim("`plugin'"))
	if "`BIAS_type'"=="" 	local BIAS_type="ds"
	if	"`BIAS_type'"!="ds" 	& 	///
		"`BIAS_type'"!="ml" 	{
		noi di in red "{bf:plugin} must be one of the following options: ds, ml"
		noi di in red "default method: ds"
		error 198
	}
	if "`BIAS_type'"=="ds" 			local COD_bias=1
	if "`BIAS_type'"=="ml" 			local COD_bias=2

	local VAR_type = lower(trim("`plugin'"))
	if "`VAR_type'"=="" 			local VAR_type="ds"
	if "`VAR_type'"=="ds" 			local COD_var=1
	if "`VAR_type'"=="ml" 			local COD_var=2
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Options for simulation-based confidence intervals
*------------------------------------------------------------------------------------------------
if "`fast'"=="" {
	if "`reps'"=="" 	local CIMC_reps=1000
	else				local CIMC_reps=`reps'
	if "`rseed'"=="" 	local CIMC_seed=-1
	else 				local CIMC_seed=`rseed'

	cap erase `"`c(tmpdir)'/STWALS_000001.tmp"'
	noi _savingopt_parse CIMC_save CIMC_savereplace : saving ".dta" `"`saving'"'
	if `"`CIMC_save'"'!="" & `"`CIMC_savereplace'"'=="" confirm new file `"`CIMC_save'"'
	if `"`CIMC_save'"'=="" {
		local CIMC_save `"`c(tmpdir)'/STWALS_000001.tmp"'
		local CIMC_savereplace "replace"
	}
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Option for variance of the error term 
*------------------------------------------------------------------------------------------------
if `sigma'<0 {
	di as err "{bf:sigma} must be nonnegative"
	error 411
}
else if `sigma'==1 {
	local mse1 "mse1"
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Display options
*------------------------------------------------------------------------------------------------
if "`fast'"=="" {
	if "`cformat'"!="" local cformat "cformat(`cformat')"
	local display_options "`table' `focus' `auxiliary' `header' `cformat'"
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Constant term (if any)
*------------------------------------------------------------------------------------------------
local nc `constant'
local ac "`auxconstant'"
if "`nc'"~="noconstant" {
	tempvar const
	gen byte `const'=1 if `touse'
}
else {
	if "`ac'"!="" {
		noi di in red "cannot specify both {bf:noconstant} and {bf:auxconstant}"
		error 184
	}
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Check for perfect collinearity in the unrestricted model
*------------------------------------------------------------------------------------------------
* collcheck (undocumented option)
if "`collcheck'"!="nocollcheck" {

	* Linear regression model 
	if `N0'>1 {
		* LS options
		local LS_opt "noheader notable `mse1'"
		
		* LS estimates of unrestricted model
		noi regress `depvar' `foclist0' `auxlist0' `wgt_exp' if `touse', `nc' `LS_opt' 
		tempname LSU_b
		matrix `LSU_b'=e(b)
		local LSU_k: colsof `LSU_b'
		local LSU_names: coln `LSU_b'

		* List of noncollinear regressors, excluding the constant term
		_ms_omit_info e(b)
		tempname x0_omit
		matrix `x0_omit'=r(omit)
		local x_names ""
		local jj=1
		foreach xx of local LSU_names {
			if `x0_omit'[1,`jj']==0 & "`xx'"!="_cons"	local x_names "`x_names' `xx'"
			local jj=`jj'+1
		}
		
		* List of noncollinear focus and auxiliary regressors, excluding the constant term 
		fvexpand `foclist0'
		local foclist0_exp `r(varlist)'
		fvexpand `auxlist0'
		local auxlist0_exp `r(varlist)'
		local a_names: list x_names-foclist0_exp  
		local f_names: list x_names-auxlist0_exp  

		* List of noncollinear focus regressors (as temp vars), including the constant term (if any)
		local MATA_foc_vars ""
		foreach xx of local f_names {
			fvrevar `xx'
			local MATA_foc_vars 	"`MATA_foc_vars' `r(varlist)'"
		}
		if "`nc'"~="noconstant" & "`ac'"=="" 	local MATA_foc_vars `const' `MATA_foc_vars' 
		local MATA_foc_names 	`f_names'
		if "`nc'"~="noconstant" & "`ac'"==""	local MATA_foc_names _cons `MATA_foc_names'

		* List of noncollinear auxiliary regressors (as temp vars), including the constant term (if any)
		local MATA_aux_vars ""
		foreach xx of local a_names {
			fvrevar `xx'
			local MATA_aux_vars 	"`MATA_aux_vars' `r(varlist)'"
		}
		if "`nc'"~="noconstant" & "`ac'"!="" 	local MATA_aux_vars  `MATA_aux_vars' `const' 
		local MATA_aux_names 	`a_names'
		if "`nc'"~="noconstant" & "`ac'"!=""	local MATA_aux_names `MATA_aux_names' _cons
		
		* Estimation sample from the unrestricted model
		replace `touse'=e(sample)
		sum `touse' `wgt_exp' if `touse' , meanonly
		local N=r(sum) 
		if `N'==0 	error 2000
	}
	* Normal location model (with 1 obs)
	else {
		sum `depvar' if `touse'
		local N=r(N)
		tempname LSU_b 
		matrix `LSU_b'=r(mean)
		local LSU_k: colsof `LSU_b'
		local LSU_names "_cons"
		local MATA_foc_vars  "" 
		local MATA_foc_names ""
		local MATA_aux_vars  `const' 
		local MATA_aux_names "_cons"
	}
}
else {

	* Size of the estimation sample
	sum `touse' `wgt_exp' if `touse' , meanonly
	local N=r(sum)
	if `N'==0 	error 2000

	* LSU_k and LSU_names
	local LSU_k: word count `foclist0' `auxlist0' 
	local LSU_names `foclist0' `auxlist0'
	if "`nc'"~="noconstant" {
		local LSU_k=`LSU_k'+1
		local LSU_names "`LSU_names' _cons" 
	}
	
	* The constant term is a focus regressor
	if "`nc'"~="noconstant" & "`ac'"=="" {
		local MATA_foc_vars `const' `foclist0' 
		local MATA_foc_names _cons 	`foclist0'
		local MATA_aux_vars  `auxlist0'  
		local MATA_aux_names `auxlist0'
	}
	* The constant term is an auxiliary regressor
	if "`nc'"~="noconstant" & "`ac'"!="" {
		local MATA_foc_vars  `foclist0' 
		local MATA_foc_names `foclist0'
		local MATA_aux_vars  `auxlist0' `const' 
		local MATA_aux_names `auxlist0' _cons
	}
	* No constant term
	if "`nc'"=="noconstant"  {
		local MATA_foc_vars  `foclist0' 
		local MATA_foc_names `foclist0'
		local MATA_aux_vars  `auxlist0'  
		local MATA_aux_names `auxlist0' 
	}
}	

* List of all noncollinear regressors (as temp vars), including the constant term (if any)
local MATA_names `MATA_foc_names' `MATA_aux_names'
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Number of focus and auxiliary regressors 
*------------------------------------------------------------------------------------------------
local k1: word count `MATA_foc_vars' 
local k2: word count `MATA_aux_vars' 
if `k2'<1 {
    noi di 
	noi di as err "The model must contain at least one auxiliary regressor (including, if any, the constant term)" 
	error 102
}
local k = `k1'+`k2'
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Check: N\ge K  
*------------------------------------------------------------------------------------------------
if `N'<`k' {
	noi di as err "Number of obs cannot be lower than the number of focus and auxiliary regressors" 
	error 198
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* WALS Estimation Procedure (in MATA)
*------------------------------------------------------------------------------------------------
noi mata: WALSMATA("`touse'","`depvar'","`MATA_foc_vars'","`MATA_aux_vars'","`prior_par'","`wgt_var'") 	
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Sample size, residual dof and standard deviation of depvar 
*------------------------------------------------------------------------------------------------
local N=int(`N')
local df_r=`N'-`k'
tempname stddev 
scalar `stddev'=r(s)[1,1]
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Estimated coefficients and their moments  
*------------------------------------------------------------------------------------------------
tempname M0_b 
matrix `M0_b'		=r(b)
if matmissing(`M0_b')==1 {
	tempname SO_Tranf_err
	scalar `SO_Tranf_err'=r(SOT_err)[1,1]
	if `SO_Tranf_err'==1 {
		noi di as err "Numerical problems in the semiortogonal transformation"
	}	
	else {
		noi di as err "WALS estimates contains missing values"
	}
	error 504
}
if "`fast'"=="" {
	tempname M0_bias M0_V M0_MSE M0_rmse  
	matrix `M0_bias'	=r(b_B)
	matrix `M0_V'		=r(b_V)
	matrix `M0_MSE'		=r(b_MSE)
	matrix `M0_rmse'	=r(b_RMSE)
}
tempname M_t_ratios M_pm M_W
matrix `M_t_ratios'=r(t_ratios) 
matrix `M_pm'=r(pm)
matrix `M_W'=r(W)
matrix coln `M_t_ratios'= `MATA_aux_names'
matrix coln `M_pm' 		= `MATA_aux_names'
matrix coln `M_W' 		= `MATA_aux_names'

* If any, put constant term as last element  
tempname M_b 
if "`nc'"~="noconstant" & "`ac'"=="" {
	local M_names: subinstr local MATA_names "_cons" ""
	local M_names `M_names' _cons
	local M_foc_names: subinstr local MATA_foc_names "_cons" ""
	local M_foc_names `M_foc_names' _cons
	local M_aux_names `MATA_aux_names'
	matrix `M_b'	=(`M0_b'[1,2..`k'],`M0_b'[1,1])
	matrix coln `M_b' 		= `M_names'
	if "`fast'"=="" {
		tempname M_bias M_V M_MSE M_rmse  
		foreach ss in bias rmse {
			matrix `M_`ss''	=(`M0_`ss''[1,2..`k'],`M0_`ss''[1,1])
			matrix coln `M_`ss'' 	= `M_names'
		}
		foreach ss in V MSE {
			matrix `M_`ss''	=(`M0_`ss''[2..`k',2..`k'],`M0_`ss''[2..`k',1])\(`M0_`ss''[1,2..`k'],`M0_`ss''[1,1])
			matrix coln `M_`ss'' 	= `M_names'
			matrix rown `M_`ss'' 	= `M_names'
		}
	}
}
else {
    local M_names `MATA_names'
	local M_foc_names `MATA_foc_names'
	local M_aux_names `MATA_aux_names'
	matrix `M_b'				=`M0_b'
	matrix coln `M_b' 			=`M_names'
	if "`fast'"=="" {
		tempname M_bias M_V M_MSE M_rmse 
		foreach ss in bias rmse {
			matrix `M_`ss''			=`M0_`ss''
			matrix coln `M_`ss'' 	=`M_names'
		}
		foreach ss in V MSE {
			matrix `M_`ss''			=`M0_`ss''
			matrix coln `M_`ss'' 	=`M_names'
			matrix rown `M_`ss'' 	=`M_names'
		}
	}
}

* Expand estimated objects to include omitted regressors 
tempname WALS_b 
matrix `WALS_b'=J(1,`LSU_k',0)
local ss=1
forvalues jj=1(1)`LSU_k'{
	local xx: word `jj' of `LSU_names'
	local yy: word `ss' of `M_names'
	if "`xx'"=="`yy'" {
		matrix `WALS_b'[1,`jj']=`M_b'[1,`ss']
		local ++ss
	}
}
matrix coln `WALS_b'=`LSU_names'
if "`fast'"=="" {
	tempname WALS_bias WALS_V WALS_MSE WALS_rmse
	foreach mm in bias rmse {
		matrix `WALS_`mm''=J(1,`LSU_k',0)
		local ss=1
		forvalues jj=1(1)`LSU_k'{
			local xx: word `jj' of `LSU_names'
			local yy: word `ss' of `M_names'
			if "`xx'"=="`yy'" {
				matrix `WALS_`mm''[1,`jj']=`M_`mm''[1,`ss']
				local ++ss
			}
		}
		matrix coln `WALS_`mm''=`LSU_names'
	}
	foreach mm in V MSE {
		matrix `WALS_`mm''=J(`LSU_k',`LSU_k',0)
		local s_rr=1
		forvalues rr=1(1)`LSU_k'{
			local xx_rr: word `rr'   of `LSU_names'
			local yy_rr: word `s_rr' of `M_names'
			local s_jj=1
			forvalues jj=1(1)`LSU_k'{
				local xx_jj: word `jj'   of `LSU_names'
				local yy_jj: word `s_jj' of `M_names'
				if "`xx_rr'"=="`yy_rr'" & "`xx_jj'"=="`yy_jj'" {
					matrix `WALS_`mm''[`rr',`jj']=`M_`mm''[`s_rr',`s_jj']
					local ++s_jj
				}
			}
			if "`xx_rr'"=="`yy_rr'" local ++s_rr
		}
		matrix coln `WALS_`mm''=`LSU_names'
		matrix rown `WALS_`mm''=`LSU_names'
	}
}

* Initialize confidence intervals
if "`fast'"=="" {
	tempname WALS_CI	
	matrix `WALS_CI'=J(2,`LSU_k',.)
	matrix rown `WALS_CI'= ll ul
	matrix coln `WALS_CI'=`LSU_names'
}

* List of omitted regressors (only base categories, not other collinear regressors)
local M_focaux_names `M_foc_names' `M_aux_names' 
local OMITTED_vars: list LSU_names - M_focaux_names
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Finalize dataset containing the MC draws of bias-corrected WALS estimates
*------------------------------------------------------------------------------------------------
if "`fast'"=="" {
	preserve
	use `"`CIMC_save'"', clear
	local jj=1
	foreach xx in `MATA_foc_names' `MATA_aux_names' {
		if "`nc'"~="noconstant" & "`ac'"=="" {
			if `jj'==1 	{
				rename v`jj' wals_bc_`k'
				lab var wals_bc_`k' "`xx'" 
			}	
			else {
				local ll=`jj'-1
				rename v`jj' wals_bc_`ll'
				lab var wals_bc_`ll' "`xx'" 
			}
		}
		else {
			rename v`jj' wals_bc_`jj'
			lab var wals_bc_`jj' "`xx'" 
		}
		local jj=`jj'+1
	}
	if "`nc'"~="noconstant" & "`ac'"=="" order wals_bc_`k', last
	save `"`CIMC_save'"', `CIMC_savereplace'
	local wals_bc_fname `"`c(filename)'"'
	restore
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Return estimation results
*------------------------------------------------------------------------------------------------
if "`fast'"!="" {
	ereturn post `WALS_b' `wgt_exp_post', dep(`depvar') obs(`N') esample(`touse') dof(`df_r') buildfvinfo 
	ereturn scalar k0 		=`LSU_k'
	ereturn scalar k1 		=`k1'
	ereturn scalar k2 		=`k2'
	ereturn scalar sigma	=`stddev'
	ereturn local title 	"WALS estimates"
	if "`PRI'"=="laplace" {
		ereturn scalar quadnpts	= .
		ereturn scalar quadatol	= .
		ereturn scalar quadrtol = .
	}
	else {
		if "`quadm'"=="gauss" {
			ereturn scalar quadnpts	=`quadnpts'
			ereturn scalar quadatol	=.
			ereturn scalar quadrtol	=.
		}
		else {
			ereturn scalar quadnpts	=.
			ereturn scalar quadatol	=`quadatol'
			ereturn scalar quadrtol	=`quadrtol'
		}
	}
	ereturn local predict 		"wals_p"
	if "`PRI'"=="laplace" 		ereturn local quadmethod	"analytic"
	else 						ereturn local quadmethod	"`quadm'"
	ereturn local prior 		"`PRI'"
	if "`nc'"=="noconstant"						ereturn local constype "`nc'"
	else if "`nc'"~="noconstant" & "`ac'"=="" 	ereturn local constype "focus"
	else 										ereturn local constype "auxiliary"
	ereturn local auxvars 		"`M_aux_names'"
	ereturn local focvars 		"`M_foc_names'"
	ereturn local omitvars	 	"`OMITTED_vars'"
	ereturn local allvars 		"`LSU_names'"
	
	ereturn local cmd 			"wals"
	ereturn matrix priorpar		=`prior_par'
	ereturn matrix wals_wgt		=`M_W'
	ereturn matrix wals_pm		=`M_pm'
	ereturn matrix t_ratios		=`M_t_ratios'
}
else{
	ereturn post `WALS_b' `WALS_V' `wgt_exp_post', dep(`depvar') obs(`N') esample(`touse') dof(`df_r') buildfvinfo
	ereturn scalar k0 			=`LSU_k'
	ereturn scalar k1 			=`k1'
	ereturn scalar k2 			=`k2'
	ereturn scalar rank 		=`k'
	ereturn scalar sigma		=`stddev'
	ereturn local title 		"WALS estimates"
	ereturn scalar level    	=`level'
	if "`PRI'"=="laplace" {
		ereturn scalar quadnpts	= .
		ereturn scalar quadatol	= .
		ereturn scalar quadrtol = .
	}
	else {
		if "`quadm'"=="gauss" {
			ereturn scalar quadnpts	=`quadnpts'
			ereturn scalar quadatol	=.
			ereturn scalar quadrtol	=.
		}
		else {
			ereturn scalar quadnpts	=.
			ereturn scalar quadatol	=`quadatol'
			ereturn scalar quadrtol	=`quadrtol'
		}
	}
	ereturn scalar reps 		=`CIMC_reps'

	ereturn matrix priorpar		=`prior_par'
	ereturn matrix wals_wgt		=`M_W'
	ereturn matrix wals_pm		=`M_pm'
	ereturn matrix t_ratios		=`M_t_ratios'
	ereturn matrix ci			=`WALS_CI'
	ereturn matrix MSE			=`WALS_MSE'
	ereturn matrix rmse			=`WALS_rmse'
	ereturn matrix bias			=`WALS_bias'

	ereturn local marginsok 	"XB default"
	ereturn local predict 		"wals_p"
	ereturn local bcsimdata 	`"`wals_bc_fname'"'	
	ereturn local plugin		"`BIAS_type'"
	if "`PRI'"=="laplace" 		ereturn local quadmethod	"analytic"
	else 						ereturn local quadmethod	"`quadm'"
	ereturn local prior 		"`PRI'"
	if "`nc'"=="noconstant"						ereturn local constype "`nc'"
	else if "`nc'"~="noconstant" & "`ac'"=="" 	ereturn local constype "focus"
	else 										ereturn local constype "auxiliary"
	ereturn local auxvars 		"`M_aux_names'"
	ereturn local focvars 		"`M_foc_names'"
	ereturn local omitvars	 	"`OMITTED_vars'"
	ereturn local allvars 		"`LSU_names'"
	ereturn local cmd 			"wals"
}
*------------------------------------------------------------------------------------------------




*------------------------------------------------------------------------------------------------
* Display estimation results 
*------------------------------------------------------------------------------------------------
if "`fast'"=="" 	noi Di_wals, level(`level') `display_options'   
*------------------------------------------------------------------------------------------------
	
}
end
*-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// Mata routines
//-----------------------------------------------------------------------------------------------------------------------	
version 14.0
mata:

//-----------------------------------------------------------------------------------------------------------------------	
// WALSMATA
//-----------------------------------------------------------------------------------------------------------------------	
void WALSMATA(touse,depvar,foc,aux,prior_par,wgtv)
{

//-----------------------------------------------------------------------------------------------------------------------	
// Declarations
//-----------------------------------------------------------------------------------------------------------------------	
	real matrix 	X1, X2, 	 			 										/*
		*/			Cwn_X1X1, Cwn_X2X2, Cwn_X1X2, ICwn_X1X1, Cwn_X2M1X2,			/*
		*/			Delta1, Delta2, Psi,											/*
		*/			Z1, D2, Z2,														/*
		*/			Cwn_Z1Z1, Cwn_Z2Z2, Cwn_Z1Z2, ICwn_Z1Z1, Cwn_Z2M1Z2,			/*
		*/			Q, D12, hat_var_b1_R,											/*
		*/			hat_var_g1_W, hat_var_g2_W, hat_cov_g1g2_W,						/*
		*/			hat_var_b1_W, hat_var_b2_W, hat_cov_b1b2_W, hat_var_b_W,		/*
		*/			hat_MSE_b_W,													/*
		*/			DRW_BCPM, DRW_b2_W, DRW_b1_R, DRW_b1_W, DRW_b_W	
		

	real colvector 	y, p_par, QP, QW, 												/*
		*/			Cwn_Z1y, Cwn_Z2y, hat_g2_U, x, hat_b1_R,						/*
		*/			POST_m, hat_delta, hat_V,										/*
		*/ 			hat_g1_W, hat_g2_W, hat_b1_W, hat_b2_W, hat_b_W,				/*
		*/			hat_bias_g1_W, hat_bias_g2_W, 									/*
		*/			hat_bias_b1_W, hat_bias_b2_W, hat_bias_b_W,						/*
		*/ 			hat_rmse_b_W													
	
	real scalar 	n, k1, k2, k, p_num, 											/*
		*/			COD_QM, COD_BIAS, COD_VAR, 										/*
		*/			NQP, AT, RT, PM_cutoff, CIMC_reps, CIMC_seed, sig_fix,			/*
		*/			Cwn_yy, s, s2
		
	string 			FAST, MC_save, comm
//-----------------------------------------------------------------------------------------------------------------------	





//-----------------------------------------------------------------------------------------------------------------------	
// Step 1: Define model and type of WALS estimations 
//-----------------------------------------------------------------------------------------------------------------------	
// Step 1.a: Load data 
	n	=strtoreal(st_local("N"))
	k1	=strtoreal(st_local("k1"))
	k2	=strtoreal(st_local("k2"))
	k	=k1+k2
	y  	=st_data(., (st_varindex(tokens(depvar))), touse)
	if (k1>0) {
	    X1 	=st_data(., (st_varindex(tokens(foc))), touse)
	}
	else {
	    X1 =J(n,k1,0)
	}
	X2 	=st_data(., (st_varindex(tokens(aux))), touse)
	wgt =st_data(., (st_varindex(tokens(wgtv))), touse)

// Step 1.b: define prior and quadrature method
	p_num	=strtoreal(st_local("COD_prior"))
	p_par 	=st_matrix(prior_par)
	COD_QM	=strtoreal(st_local("COD_quadm"))

// Step 1.c: other WALS options

	// points and weights for Gauss-Laguerre/Gauss-Legendre quadrature 
	if (COD_QM==1) {
		if (st_local("quadext")!="") {
			//external GL_PW
			GL_PW=valofexternal(st_local("quadext"))
			NQP		=rows(GL_PW)
		}
		else {
			NQP		=strtoreal(st_local("quadnpts"))
			if (p_num==1)	GL_PW 	=J(2,2,0)
			else {
			    if (p_num<7)	GL_PW 	=gausslaguerre(NQP)
				else			GL_PW 	=gausslegendre(NQP)
			}
		}
		QP		=GL_PW[.,1]
		QW 		=GL_PW[.,2]
	}
	// absolute and relative tolerances for adaptive quadrature
	else {
		AT		=strtoreal(st_local("quadatol"))
		RT		=strtoreal(st_local("quadrtol"))
	}

	// cutoff point to rescale of the posterior mean
	PM_cutoff	=8.5

	// fast (point estimates only)
	FAST=st_local("fast")

	if (FAST=="") {
	    
		// plug-in estimators of bias and variance of the posterior mean
		COD_BIAS	=strtoreal(st_local("COD_bias"))
		COD_VAR		=strtoreal(st_local("COD_var"))
		
		// simulation-based confidence intervals
		CIMC_reps  	=strtoreal(st_local("CIMC_reps"))
		CIMC_seed  	=strtoreal(st_local("CIMC_seed"))
		CIMC_save  	=`"""'+st_local("CIMC_save")+`"""'
		CIMC_saverep=st_local("CIMC_savereplace")
		if (CIMC_seed>0) rseed(CIMC_seed)

		// tabulated values of bias and variance of the posterior mean
		external tabmom
	}

	// fixed standard deviation of the error term
	sig_fix	=strtoreal(st_local("sigma"))	
//-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// Step 2: scaling
//-----------------------------------------------------------------------------------------------------------------------	
	Cwn_X1X1	=makesymmetric(quadcross(X1,wgt,X1):/n) 
	Cwn_X2X2	=makesymmetric(quadcross(X2,wgt,X2):/n) 
	Cwn_X1X2	=quadcross(X1,wgt,X2):/n 
	ICwn_X1X1	=invsym(Cwn_X1X1)
	Cwn_X2M1X2	=makesymmetric(Cwn_X2X2 - Cwn_X1X2' * ICwn_X1X1 * Cwn_X1X2)
	Delta1		=diag(diag(Cwn_X1X1):^(-.5))
	Delta2		=diag(diag(Cwn_X2M1X2):^(-.5))
	Psi			=makesymmetric(Delta2*Cwn_X2M1X2*Delta2)
//-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// Step 3: Semi-orthogonalization transformation 
//-----------------------------------------------------------------------------------------------------------------------	
	Z1			=X1 * Delta1
	D2			=Delta2*matpowersym(Psi,-.5) 	
	Z2			=X2 * D2
	if (missing(D2)>0) {
		st_matrix("r(SOT_err)"			, 1	)	
	}
	else {
		st_matrix("r(SOT_err)"			, 0	)	
	}
		
	//matpowersym(Psi,-.5) *sqrt(n)	
	Cwn_Z1Z1  	=quadcross(Z1,wgt,Z1):/n
    Cwn_Z2Z2  	=quadcross(Z2,wgt,Z2):/n
    Cwn_Z1Z2  	=quadcross(Z1,wgt,Z2):/n
	ICwn_Z1Z1 	=invsym(Cwn_Z1Z1)
	Cwn_Z2M1Z2	=makesymmetric(Cwn_Z2Z2 - Cwn_Z1Z2' * ICwn_Z1Z1 * Cwn_Z1Z2)
	Q       	= ICwn_Z1Z1 * Cwn_Z1Z2
	D12 		= Delta1 * Q
    Cwn_Z1y  	=quadcross(Z1,wgt,y):/n
    Cwn_Z2y  	=quadcross(Z2,wgt,y):/n
    Cwn_yy  	=quadcross(y,wgt,y):/n
//-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// Step 4: LS estimates of (transformed) unrestricted and restricted models
//-----------------------------------------------------------------------------------------------------------------------	
// Step 4.a: unrestricted model 
	hat_g2_U		=Cwn_Z2y-Cwn_Z1Z2' * ICwn_Z1Z1 * Cwn_Z1y
	if (sig_fix==0) {
		s2			=((Cwn_yy-Cwn_Z1y'*ICwn_Z1Z1*Cwn_Z1y)-hat_g2_U'*hat_g2_U)*(n/(n-k))
		s           =sqrt(s2)
	}
	else {
		s 			=sig_fix
		s2 			=s^2
	}	
	x           	=n^.5 * hat_g2_U / s 

// Step 4.b: restricted model
	hat_b1_R		=Delta1 * ICwn_Z1Z1 * Cwn_Z1y
	hat_var_b1_R	=makesymmetric((s2/n) * Delta1 *  ICwn_Z1Z1 * Delta1 )
//-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// Step 5: posterior means in the NLM
//-----------------------------------------------------------------------------------------------------------------------	
	if (COD_QM==1)		POST_m=NLM_PM_GL(x,p_num,p_par,NQP,QP,QW,PM_cutoff)
	else				POST_m=NLM_PM_AQ(x,p_num,p_par,AT,RT,PM_cutoff)	
//-----------------------------------------------------------------------------------------------------------------------	




//-----------------------------------------------------------------------------------------------------------------------	
// Step 6: Plug-in estimators of the bias of the posterior mean 
//-----------------------------------------------------------------------------------------------------------------------	
if (FAST=="") {
	// Step 6.a: double-shrinkage 
	if (COD_BIAS==1) 	hat_delta	= PLUG_IN(POST_m, tabmom, 1, p_num, p_par) 

	// Step 6.b: maximum likelihood 
	if (COD_BIAS==2) 	hat_delta	= PLUG_IN(x		, tabmom, 1, p_num, p_par)
}
//-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// Step 7: Plug-in estimators of the variance of the posterior mean 
//-----------------------------------------------------------------------------------------------------------------------	
if (FAST=="") {
	// Step 7.a: double-shrinkage 
	if (COD_VAR==1) 	hat_V	= PLUG_IN(POST_m, tabmom, 2, p_num, p_par) 

	// Step 7.b: maximum likelihood 
	if (COD_VAR==2) 	hat_V	= PLUG_IN(x     , tabmom, 2, p_num, p_par) 
}
//-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// Step 8: WALS point estimates 
//-----------------------------------------------------------------------------------------------------------------------	
// Step 8.a: estimates of \gamma_2
hat_g2_W 		=(s/n^.5) * POST_m

// Step 8.b: estimates of \gamma_1
hat_g1_W		=ICwn_Z1Z1 * (Cwn_Z1y - Cwn_Z1Z2 * hat_g2_W)			

// Step 8.c: estimates of \beta_2
hat_b2_W		=D2 * hat_g2_W

// Step 8.d: estimates of \beta_1
hat_b1_W 		=Delta1 * hat_g1_W

// Step 8.e: estimates of \beta
hat_b_W			=hat_b1_W\hat_b2_W
//-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// Step 9: Estimated biases of the WALS estimator
//-----------------------------------------------------------------------------------------------------------------------	
if (FAST=="") {
	// Step 9.a: estimates of bias(\hat \gamma_2)
	hat_bias_g2_W 	=(s/n^.5) * hat_delta

	// Step 9.b: estimates of bias(\hat \gamma_1)
	hat_bias_g1_W	=-Q * hat_bias_g2_W			

	// Step 9.c: estimates of bias(\hat \beta_2)
	hat_bias_b2_W	=D2 * hat_bias_g2_W

	// Step 9.d: estimates of bias(\hat \beta1)
	hat_bias_b1_W	=Delta1 * hat_bias_g1_W

	// Step 9.e: estimates of bias(\hat \beta)
	hat_bias_b_W	=hat_bias_b1_W\hat_bias_b2_W
}
//-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// Step 10: Estimated variance of the WALS estimator 
//-----------------------------------------------------------------------------------------------------------------------	
if (FAST=="") {
	// Step 10.a: estimates of var(\hat \gamma_2)
	hat_var_g2_W	=(s2/n) * diag(hat_V)

	// Step 10.b: estimates of var(\hat \gamma_1)
	hat_var_g1_W  	=(s2/n) * ICwn_Z1Z1 + Q * hat_var_g2_W * Q'

	// Step 10.c: estimates of cov(\hat \gamma_1,\hat \beta_2)
	hat_cov_g1g2_W  =-Q * hat_var_g2_W

	// Step 10.d: estimates of var(\hat \beta_2)
	hat_var_b2_W	=D2 * hat_var_g2_W * D2'

	// Step 10.e: estimates of var(\hat \beta_1)
	hat_var_b1_W 	=Delta1 * hat_var_g1_W * Delta1'

	// Step 10.f: estimates of cov(\hat \beta_1,\hat \beta_2)
	hat_cov_b1b2_W  =Delta1 * hat_cov_g1g2_W * D2'

	// Step 10.g: estimates of var(\hat \beta)
	hat_var_b_W		=makesymmetric((hat_var_b1_W, hat_cov_b1b2_W)\(hat_cov_b1b2_W',hat_var_b2_W))
}
//-----------------------------------------------------------------------------------------------------------------------	


	
//-----------------------------------------------------------------------------------------------------------------------	
// Step 11: Estimated MSE matrix of the WALS estimator
//-----------------------------------------------------------------------------------------------------------------------	
if (FAST=="") {
	// Step 11.a: estimates of MSE(\hat \beta)
	hat_MSE_b_W		=hat_bias_b_W*hat_bias_b_W'+hat_var_b_W

	// Step 11.b: estimates of RMSE(\hat \beta_j), j=1...k
	hat_rmse_b_W	=(diagonal(hat_MSE_b_W)):^.5
}
//-----------------------------------------------------------------------------------------------------------------------	


		
//---------------------------------------------------------------------------------------------	
// Step 12: Draws from the (estimated) sampling distribution of the (bias-corrected) WALS estimator
//---------------------------------------------------------------------------------------------	
if (FAST=="") {
	// Step 12.a: Draws of bias-corrected estimates of \beta_2	
	if (COD_QM==1) 	{
		DRW_BCPM=DRAWS_BCPM_GL(POST_m-hat_delta, CIMC_reps, p_num, p_par, NQP, QP, QW, PM_cutoff, COD_BIAS, tabmom)
	}
	else {
		DRW_BCPM=DRAWS_BCPM_AQ(POST_m-hat_delta, CIMC_reps, p_num, p_par, AT, RT,      PM_cutoff, COD_BIAS, tabmom)
	}
	DRW_b2_W=(s/n^.5)*DRW_BCPM*D2'

	// Step 12.b: Draws of bias-corrected estimates of \beta_1 
	if (k1>0) {
		
		DRW_U1			=J(CIMC_reps,k1,0)
		primes			=FirstPrimes(k1)
		for (h=1; h<=k1; h++) {
			DRW_U1[.,h]	=ghalton(CIMC_reps,primes[1,h],runiform(1,1))
		}
		DRW_b1_R		=hat_b1_R' # J(CIMC_reps,1,1)+invnormal(DRW_U1)*cholesky(hat_var_b1_R)'
		//DRW_b1_R
		//DRW_b1_R		=hat_b1_R' # J(CIMC_reps,1,1)+rnormal(CIMC_reps,k1,0,1)*cholesky(hat_var_b1_R)'
		DRW_b1_W		=DRW_b1_R - (s/n^.5)*DRW_BCPM*D12'
	}
	else DRW_b1_W		=J(CIMC_reps,0,0)

	// Step 12.c: draws of bias-corrected estimates of \beta 
	DRW_b_W=DRW_b1_W,DRW_b2_W
		
	// Step 12.d: Store draws of bias-corrected estimates 
	comm	="preserve"
	stata(comm, 1)
	st_dropvar(.)
	for (i=1; i<=k; i++) {
		(void) st_addvar("double", sprintf("v%g",i))
	}
	st_addobs(CIMC_reps)
	st_store(.,., DRW_b_W)
	comm	="save " + CIMC_save + " , " + CIMC_saverep 
	stata(comm, 1)
	comm	="restore"
	stata(comm, 1)
}
//---------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// Step 13: Return results  
//-----------------------------------------------------------------------------------------------------------------------	
st_matrix("r(b)"			, hat_b_W'			)	
st_matrix("r(s)"			, s					)
st_matrix("r(t_ratios)"		, x'				)	
st_matrix("r(pm)"			, POST_m'			)	
st_matrix("r(W)"			, (POST_m:/x)'		)	
if (FAST=="") {
	st_matrix("r(b_B)"		, hat_bias_b_W'		)	
	st_matrix("r(b_V)"		, hat_var_b_W		)	
	st_matrix("r(b_MSE)"	, hat_MSE_b_W	 	)
	st_matrix("r(b_RMSE)"	, hat_rmse_b_W' 	)
}
//-----------------------------------------------------------------------------------------------------------------------	

}		
//-----------------------------------------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------------------------------	
// Routine			: I_INT_AQ1 
// Description		: Integral I(x) by adaptive quadrature 
//-----------------------------------------------------------------------------------------------------------------------	
function I_INT_f1(real scalar t, real scalar x) {
	return(exp(x*(1-t))/t) 
}
real matrix I_INT_AQ1(real colvector x)
{

	// Number of data points to be evaluated
	dim_x=rows(x)

	// Define integration problem 
	class Quadrature scalar   I_INT1
	I_INT1 = Quadrature()
	I_INT1.setEvaluator(&I_INT_f1())
	I_INT1.setLimits((1, .))
	//I_INT1.setAbstol(1e-13)
	//I_INT1.setReltol(1e-11)
	//I_INT1.setAbstol(1e-7)
	//I_INT1.setReltol(1e-5)

	// Loop over values of x & integration
	RES=J(dim_x,1,0)
	for (nn=1; nn<=dim_x; nn++) {
		I_INT1.setArgument(1, x[nn])
		RES[nn,1]=I_INT1.integrate()
	}
	
	// Return results
	return(RES) 
}
//-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// Routine			: NLM_PRIOR
// Description		: Prior densities for the normal location model
//-----------------------------------------------------------------------------------------------------------------------	
function NLM_PRIOR(real colvector theta, real scalar p_num,|real colvector p_par) {

	// Laplace
	if (p_num==1) {
		pp_c=p_par[1]
		return(laplaceden(0,1/pp_c,theta))
	}
	// Subbotin
	if (p_num==2) {
		pp_b=p_par[1]
		pp_c=p_par[2]
		return(((pp_b :* pp_c:^(1/pp_b) :* exp(-lngamma(1/pp_b))):/2) :* exp(-pp_c:*(abs(theta):^(pp_b))))
	}
	// Weibull
	if (p_num==3) {
		pp_b=p_par[1]
		pp_c=p_par[2]
		return(((pp_b :* pp_c):/2) :* abs(theta):^(pp_b-1) :* exp(-pp_c:*(abs(theta):^(pp_b))))
	}
	// Pareto
	if (p_num==4) {
		pp_a=p_par[1]
		pp_c=p_par[2]
		return((pp_c * (1-pp_a) * (1:+pp_c*abs(theta)):^(-1/pp_a)):/(2*pp_a))
	}
	// Cauchy
	if (p_num==5) {
		return(cauchyden(0,1,theta))
	}
	// Log
	if (p_num==6) {
	    pp_c=p_par[1]
		return((log(1:+(pp_c:/theta):^2)):/(2*pi()*pp_c))
	}
	// Horseshoe
	if (p_num==7) {
	    pp_c=p_par[1]
		s=.5*(theta:/pp_c):^2
		return(I_INT_AQ1(s):/(pp_c:*(2*pi()^3)^(.5)))
	}
}
//-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// Routine			: NLM_PM_GL
// Description		: Posterior mean by Gauss-Laguerre/Gauss-Legendre quadrature with rescaling for x>cutoff
//-----------------------------------------------------------------------------------------------------------------------	
real matrix NLM_PM_GL(real colvector x, real scalar p_num, real colvector p_par, real scalar NQP, real colvector QP, real colvector QW, real scalar cutoff)
{
	// Laplace
	if (p_num==1) {	
		pp_c=p_par[1]
		y=abs(x)
		psi_y=normal(-y:-pp_c):/normal(y:-pp_c)
		h_y= (1:-exp(2:*pp_c:*y):*psi_y):/(1:+exp(2:*pp_c:*y):*psi_y)
		_editmissing(h_y,1)
		PM=(x!=0):*(x:-(sign(x)):*pp_c:* h_y)
	}
	
	// Weibull, Subbotin, Pareto, Cauchy and log-type 
	if (p_num>=2 & p_num<=6) {	
		
		// Number of data points to be evaluated
		dim_x=rows(x)
		PM=J(dim_x,1,.)
			
		// Prior (only for points below the cutoff)
		if (colsum(abs(x):<=cutoff)>=1) {
			QP_Prior=NLM_PRIOR(QP, p_num, p_par) 
		}
			
		// Subbotin
		if (p_num==2) {
			pp_a=0
			pp_b=p_par[1]
			pp_c=p_par[2]
			for (nn=1; nn<=dim_x; nn++) {
				if (x[nn]==0) PM[nn,.]=0
				else {
					A=J(NQP,2,0)
					if (abs(x[nn])<=cutoff) {
						A[.,1]=(               normalden(x[nn]:-QP) +                   normalden(x[nn]:+QP)) :* QP_Prior :* exp(QP)
						A[.,2]=((x[nn]:-QP) :* normalden(x[nn]:-QP)	+ (x[nn]:+QP)	:* 	normalden(x[nn]:+QP)) :* QP_Prior :* exp(QP)
					}
					else {
						Q_x =1:+(QP:/x[nn])
						Q_mx=1:-(QP:/x[nn])
						G_x =((abs(Q_x) ):^(-pp_a)) :* (exp(-pp_c:*(abs(x[nn])):^pp_b:*((abs(Q_x) ):^pp_b:-1)))
						G_mx=((abs(Q_mx)):^(-pp_a)) :* (exp(-pp_c:*(abs(x[nn])):^pp_b:*((abs(Q_mx)):^pp_b:-1)))
						A[.,1]=     (G_mx:+G_x) :* normalden(QP) :* exp(QP)
						A[.,2]=QP:* (G_mx:-G_x) :* normalden(QP) :* exp(QP)
					}
					A=quadcolsum(A:*QW)
					if (A[1]>0 & abs(A[2]/A[1])<abs(x[nn])) PM[nn,.]=x[nn] - sign(A[2])*max((abs(A[2]/A[1]),0))
				}
			}
		}
		
		// Weibull
		if (p_num==3) {
			pp_b=p_par[1]
			pp_c=p_par[2]
			pp_a=1-pp_b
			for (nn=1; nn<=dim_x; nn++) {
				if (x[nn]==0) PM[nn,.]=0
				else {
					A=J(NQP,2,0)
					if (abs(x[nn])<=cutoff) {
						A[.,1]=(               normalden(x[nn]:-QP) +                   normalden(x[nn]:+QP)) :* QP_Prior :* exp(QP)
						A[.,2]=((x[nn]:-QP) :* normalden(x[nn]:-QP)	+ (x[nn]:+QP)	:* 	normalden(x[nn]:+QP)) :* QP_Prior :* exp(QP)
					}
					else {
						Q_x =1:+(QP:/x[nn])
						Q_mx=1:-(QP:/x[nn])
						G_x =((abs(Q_x) ):^(-pp_a)) :* (exp(-pp_c:*(abs(x[nn])):^pp_b:*((abs(Q_x) ):^pp_b:-1)))
						G_mx=((abs(Q_mx)):^(-pp_a)) :* (exp(-pp_c:*(abs(x[nn])):^pp_b:*((abs(Q_mx)):^pp_b:-1)))
						A[.,1]=     (G_mx:+G_x) :* normalden(QP) :* exp(QP)
						A[.,2]=QP:* (G_mx:-G_x) :* normalden(QP) :* exp(QP)
					}
					A=quadcolsum(A:*QW)
					if (A[1]>0 & abs(A[2]/A[1])<abs(x[nn])) PM[nn,.]=x[nn] - sign(A[2])*max((abs(A[2]/A[1]),0))
				}
			}
		}
		
		// Pareto
		if (p_num==4) {
			pp_a=p_par[1]
			pp_c=p_par[2]
			for (nn=1; nn<=dim_x; nn++) {
				if (x[nn]==0) PM[nn,.]=0
				else {
					A=J(NQP,2,0)
					if (abs(x[nn])<=cutoff) {
						A[.,1]=(               normalden(x[nn]:-QP) +                   normalden(x[nn]:+QP)) :* QP_Prior :* exp(QP)
						A[.,2]=((x[nn]:-QP) :* normalden(x[nn]:-QP)	+ (x[nn]:+QP)	:* 	normalden(x[nn]:+QP)) :* QP_Prior :* exp(QP)
					}
					else {
						G_x =((1:+pp_c:*abs(x[nn]:+QP)):/(1+pp_c*abs(x[nn]))):^(-1/pp_a)	
						G_mx=((1:+pp_c:*abs(x[nn]:-QP)):/(1+pp_c*abs(x[nn]))):^(-1/pp_a)	
						A[.,1]=     (G_mx:+G_x) :* normalden(QP) :* exp(QP)
						A[.,2]=QP:* (G_mx:-G_x) :* normalden(QP) :* exp(QP)
					}
					A=quadcolsum(A:*QW)
					if (A[1]>0 & abs(A[2]/A[1])<abs(x[nn])) PM[nn,.]=x[nn] - sign(A[2])*max((abs(A[2]/A[1]),0))
				}
			}
		}
		
		// Cauchy
		if (p_num==5) {
			for (nn=1; nn<=dim_x; nn++) {
				if (x[nn]==0) PM[nn,.]=0
				else {
					A=J(NQP,2,0)
					if (abs(x[nn])<=cutoff) {
						A[.,1]=(               normalden(x[nn]:-QP) +                   normalden(x[nn]:+QP)) :* QP_Prior :* exp(QP)
						A[.,2]=((x[nn]:-QP) :* normalden(x[nn]:-QP)	+ (x[nn]:+QP)	:* 	normalden(x[nn]:+QP)) :* QP_Prior :* exp(QP)
					}
					else {
						G_x =(1+(x[nn])^2):/(1:+(x[nn]:+QP):^2)
						G_mx=(1+(x[nn])^2):/(1:+(x[nn]:-QP):^2)
						A[.,1]=     (G_mx:+G_x) :* normalden(QP) :* exp(QP)
						A[.,2]=QP:* (G_mx:-G_x) :* normalden(QP) :* exp(QP)
					}
					A=quadcolsum(A:*QW)
					if (A[1]>0 & abs(A[2]/A[1])<abs(x[nn])) PM[nn,.]=x[nn] - sign(A[2])*max((abs(A[2]/A[1]),0))
				}
			}
		}
		
		// Log
		if (p_num==6) {
			pp_c=p_par[1]
			for (nn=1; nn<=dim_x; nn++) {
				if (x[nn]==0) PM[nn,.]=0
				else {
					A=J(NQP,2,0)
					if (abs(x[nn])<=cutoff) {
						A[.,1]=(               normalden(x[nn]:-QP) +                   normalden(x[nn]:+QP)) :* QP_Prior :* exp(QP)
						A[.,2]=((x[nn]:-QP) :* normalden(x[nn]:-QP)	+ (x[nn]:+QP)	:* 	normalden(x[nn]:+QP)) :* QP_Prior :* exp(QP)
					}
					else {
						G_x =(log(1:+pp_c^2:/(x[nn]:+QP):^2)):/(log(1:+pp_c^2/x[nn]^2))	
						G_mx=(log(1:+pp_c^2:/(x[nn]:-QP):^2)):/(log(1:+pp_c^2/x[nn]^2))		
						A[.,1]=     (G_mx:+G_x) :* normalden(QP) :* exp(QP)
						A[.,2]=QP:* (G_mx:-G_x) :* normalden(QP) :* exp(QP)
					}
					A=quadcolsum(A:*QW)
					if (A[1]>0 & abs(A[2]/A[1])<abs(x[nn])) PM[nn,.]=x[nn] - sign(A[2])*max((abs(A[2]/A[1]),0))
				}
			}
		}
	}

	// Horseshoe
	if (p_num==7) {
		dim_x=rows(x)
		PM=J(dim_x,1,.)
		u=.5:*(QP:+1)
		pp_c=p_par[1]
		for (nn=1; nn<=dim_x; nn++) {
			if (x[nn]==0|abs(x[nn])>=600) PM[nn,.]=x[nn]
			else {
				A=J(NQP,2,0)
				if (abs(x[nn])<=cutoff) {
					A[.,1]=   (1:-u):^(-.5):*exp(-x[nn]^2:*u:/2):/(1:+(pp_c^2-1):*u) 
					A[.,2]=u:*(1:-u):^(-.5):*exp(-x[nn]^2:*u:/2):/(1:+(pp_c^2-1):*u) 
				}
				else {
					A[.,1]=x[nn]^2:*   (1:-u):^(-.5):*exp(-x[nn]^2:*u:/2):/(1:+(pp_c^2-1):*u) 
					A[.,2]=x[nn]^2:*u:*(1:-u):^(-.5):*exp(-x[nn]^2:*u:/2):/(1:+(pp_c^2-1):*u) 
				}
				A=.5:*quadcolsum(A:*QW)
				if (A[2]/A[1]>0 & A[2]/A[1]<1) PM[nn,.]=(1-A[2]/A[1])*x[nn]
			}
		}
	}

	// Return results on the posterior mean
	return(PM) 
}
//-----------------------------------------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------------------------------	
// Routine			: NLM_PM_AQ
// Description		: Posterior mean by adaptive quadrature with rescaling for x>cutoff
//-----------------------------------------------------------------------------------------------------------------------	
// Subbotin evaluator
function PM_AQ_A0_2_f(real scalar theta, real scalar x, real colvector p_par, real scalar cutoff) {
    if (abs(x)<=cutoff|abs(theta)==abs(x)) {
		return((          normalden(x-theta)+          normalden(x+theta))*NLM_PRIOR(theta,2,p_par)) 
	}
	else {
		pp_a=0
		pp_b=p_par[1]
		pp_c=p_par[2]
		Q_x =1+theta/x
		Q_mx=1-theta/x
		G_x =(abs(Q_x )^(-pp_a)) * (exp(-pp_c*abs(x)^(pp_b)*(abs(Q_x )^(pp_b)-1)))
		G_mx=(abs(Q_mx)^(-pp_a)) * (exp(-pp_c*abs(x)^(pp_b)*(abs(Q_mx)^(pp_b)-1)))
		return(       (G_mx+G_x) * normalden(theta))
	}
}
function PM_AQ_A1_2_f(real scalar theta, real scalar x, real colvector p_par, real scalar cutoff) {
    if (abs(x)<=cutoff|abs(theta)==abs(x)) {
		return(((x-theta)*normalden(x-theta)+(x+theta)*normalden(x+theta))*NLM_PRIOR(theta,2,p_par)) 
	}
	else {
		pp_a=0
		pp_b=p_par[1]
		pp_c=p_par[2]
		Q_x =1+theta/x
		Q_mx=1-theta/x
		G_x =(abs(Q_x )^(-pp_a)) * (exp(-pp_c*abs(x)^(pp_b)*(abs(Q_x )^(pp_b)-1)))
		G_mx=(abs(Q_mx)^(-pp_a)) * (exp(-pp_c*abs(x)^(pp_b)*(abs(Q_mx)^(pp_b)-1)))
		return(theta*(G_mx-G_x)  * normalden(theta))
	}
}
// Weibull evaluator
function PM_AQ_A0_3_f(real scalar theta, real scalar x, real colvector p_par, real scalar cutoff) {
    if (abs(x)<=cutoff|abs(theta)==abs(x)) {
		return((          normalden(x-theta)+          normalden(x+theta))*NLM_PRIOR(theta,2,p_par)) 
	}
	else {
		pp_b=p_par[1]
		pp_c=p_par[2]
		pp_a=1-pp_b
		Q_x =1+theta/x
		Q_mx=1-theta/x
		G_x =(abs(Q_x )^(-pp_a)) * (exp(-pp_c*abs(x)^(pp_b)*(abs(Q_x )^(pp_b)-1)))
		G_mx=(abs(Q_mx)^(-pp_a)) * (exp(-pp_c*abs(x)^(pp_b)*(abs(Q_mx)^(pp_b)-1)))
		return(       (G_mx+G_x) * normalden(theta))
	}
}
function PM_AQ_A1_3_f(real scalar theta, real scalar x, real colvector p_par, real scalar cutoff) {
    if (abs(x)<=cutoff|abs(theta)==abs(x)) {
		return(((x-theta)*normalden(x-theta)+(x+theta)*normalden(x+theta))*NLM_PRIOR(theta,2,p_par)) 
	}
	else {
		pp_b=p_par[1]
		pp_c=p_par[2]
		pp_a=1-pp_b
		Q_x =1+theta/x
		Q_mx=1-theta/x
		G_x =(abs(Q_x )^(-pp_a)) * (exp(-pp_c*abs(x)^(pp_b)*(abs(Q_x )^(pp_b)-1)))
		G_mx=(abs(Q_mx)^(-pp_a)) * (exp(-pp_c*abs(x)^(pp_b)*(abs(Q_mx)^(pp_b)-1)))
		return(theta*(G_mx-G_x)  * normalden(theta))
	}
}
// Pareto evaluator
function PM_AQ_A0_4_f(real scalar theta, real scalar x, real colvector p_par, real scalar cutoff) {
    if (abs(x)<=cutoff) {
		return((          normalden(x-theta)+          normalden(x+theta))*NLM_PRIOR(theta,4,p_par)) 
	}
	else {
		pp_a=p_par[1]
		pp_c=p_par[2]
		G_x =((1+pp_c*abs(x+theta))/(1+pp_c*abs(x)))^(-1/pp_a)	
		G_mx=((1+pp_c*abs(x-theta))/(1+pp_c*abs(x)))^(-1/pp_a)	
		return(       (G_mx+G_x) * normalden(theta))
	}
}
function PM_AQ_A1_4_f(real scalar theta, real scalar x, real colvector p_par, real scalar cutoff) {
    if (abs(x)<=cutoff) {
		return(((x-theta)*normalden(x-theta)+(x+theta)*normalden(x+theta))*NLM_PRIOR(theta,4,p_par)) 
	}
	else {
		pp_a=p_par[1]
		pp_c=p_par[2]
		G_x =((1+pp_c*abs(x+theta))/(1+pp_c*abs(x)))^(-1/pp_a)	
		G_mx=((1+pp_c*abs(x-theta))/(1+pp_c*abs(x)))^(-1/pp_a)	
		return(theta*(G_mx-G_x) * normalden(theta))
	}
}
// Cauchy evaluator
function PM_AQ_A0_5_f(real scalar theta, real scalar x, real scalar cutoff) {
    if (abs(x)<=cutoff) {
		return((          normalden(x-theta)+          normalden(x+theta))*NLM_PRIOR(theta,5)) 
	}
	else {
		G_x =(1+x^2)/(1+(x+theta)^2)
		G_mx=(1+x^2)/(1+(x-theta)^2)
		return(       (G_mx+G_x) * normalden(theta))
	}
}
function PM_AQ_A1_5_f(real scalar theta, real scalar x, real scalar cutoff) {
    if (abs(x)<=cutoff) {
		return(((x-theta)*normalden(x-theta)+(x+theta)*normalden(x+theta))*NLM_PRIOR(theta,5)) 
	}
	else {
		G_x =(1+x^2)/(1+(x+theta)^2)
		G_mx=(1+x^2)/(1+(x-theta)^2)
		return(theta*(G_mx-G_x) * normalden(theta))
	}
}
// Log evaluator
function PM_AQ_A0_6_f(real scalar theta, real scalar x, real scalar pp_c, real scalar cutoff) {
    if (abs(x)<=cutoff) {
		return((          normalden(x-theta)+          normalden(x+theta))*NLM_PRIOR(theta,6,pp_c)) 
	}
	else {
		G_x =(log(1+pp_c^2/(x+theta)^2))/(log(1+pp_c^2/x^2))		
		G_mx=(log(1+pp_c^2/(x-theta)^2))/(log(1+pp_c^2/x^2))
		return(       (G_mx+G_x) * normalden(theta))
	}
}
function PM_AQ_A1_6_f(real scalar theta, real scalar x, real scalar pp_c, real scalar cutoff) {
    if (abs(x)<=cutoff) {
		return(((x-theta)*normalden(x-theta)+(x+theta)*normalden(x+theta))*NLM_PRIOR(theta,6,pp_c)) 
	}
	else {
		G_x =(log(1+pp_c^2/(x+theta)^2))/(log(1+pp_c^2/x^2))		
		G_mx=(log(1+pp_c^2/(x-theta)^2))/(log(1+pp_c^2/x^2))
		return(theta*(G_mx-G_x) * normalden(theta))
	}
}
// Horseshoe evaluator
function PM_AQ_A0_7_f(real scalar u, real scalar x, real scalar pp_c, real scalar cutoff) {
    if (abs(x)<=cutoff) {
		return(      (1-u)^(-.5)*exp(-x^2*u/2)/(1+(pp_c^2-1)*u)) 
	}
	else {
		return(x^2*  (1-u)^(-.5)*exp(-x^2*u/2)/(1+(pp_c^2-1)*u)) 
	}
}
function PM_AQ_A1_7_f(real scalar u, real scalar x, real scalar pp_c, real scalar cutoff) {
    if (abs(x)<=cutoff) {
		return(    u*(1-u)^(-.5)*exp(-x^2*u/2)/(1+(pp_c^2-1)*u)) 
	}
	else {
		return(x^2*u*(1-u)^(-.5)*exp(-x^2*u/2)/(1+(pp_c^2-1)*u)) 
	}
}

// Main routine for adaptive quadrature
real matrix NLM_PM_AQ(real colvector x, real scalar p_num, real colvector p_par, real scalar AT, real scalar RT, real scalar cutoff)
{
	// Laplace
	if (p_num==1) {	
		pp_c=p_par[1]
		y=abs(x)
		psi_y=normal(-y:-pp_c):/normal(y:-pp_c)
		h_y= (1:-exp(2:*pp_c:*y):*psi_y):/(1:+exp(2:*pp_c:*y):*psi_y)
		_editmissing(h_y,1)
		PM=(x!=0):*(x:-(sign(x)):*pp_c:* h_y)
	}
	
	// Other priors
	if (p_num>=2 & p_num<=7) {	
	    
		dim_x=rows(x)
		PM=J(dim_x,1,.)
		A=J(1,2,0)

		class Quadrature scalar   PM_AQ_A0
		PM_AQ_A0 = Quadrature()
		PM_AQ_A0.setAbstol(AT)
		PM_AQ_A0.setReltol(RT)

		class Quadrature scalar   PM_AQ_A1
		PM_AQ_A1 = Quadrature()
		PM_AQ_A1.setAbstol(AT)
		PM_AQ_A1.setReltol(RT)
				
		// Subbotin
		if (p_num==2) {
			PM_AQ_A0.setEvaluator(&PM_AQ_A0_2_f())
			PM_AQ_A0.setArgument(2, p_par)
			PM_AQ_A0.setArgument(3, cutoff)

			PM_AQ_A1.setEvaluator(&PM_AQ_A1_2_f())
			PM_AQ_A1.setArgument(2, p_par)
			PM_AQ_A1.setArgument(3, cutoff)
			
			PM_AQ_A0.setLimits((0, .))
			PM_AQ_A1.setLimits((0, .))
			for (nn=1; nn<=dim_x; nn++) {
				if (x[nn]==0) PM[nn,.]=0
				else {
					PM_AQ_A0.setArgument(1, x[nn])
					A[1,1]=PM_AQ_A0.integrate()
					
					PM_AQ_A1.setArgument(1, x[nn])
					A[1,2]=PM_AQ_A1.integrate()
					
					if (PM_AQ_A0.errorcode()==0 & PM_AQ_A1.errorcode()==0 & A[1]>0 & abs(A[2]/A[1])<abs(x[nn])) {
						PM[nn,.]=x[nn] - sign(A[2])*max((abs(A[2]/A[1]),0))
					}
				}
			}
		}

		// Weibull
		if (p_num==3) {
			PM_AQ_A0.setEvaluator(&PM_AQ_A0_3_f())
			PM_AQ_A0.setArgument(2, p_par)
			PM_AQ_A0.setArgument(3, cutoff)

			PM_AQ_A1.setEvaluator(&PM_AQ_A1_3_f())
			PM_AQ_A1.setArgument(2, p_par)
			PM_AQ_A1.setArgument(3, cutoff)
			
			PM_AQ_A0.setLimits((0, .))
			PM_AQ_A1.setLimits((0, .))
			for (nn=1; nn<=dim_x; nn++) {
				if (x[nn]==0) PM[nn,.]=0
				else {
					PM_AQ_A0.setArgument(1, x[nn])
					A[1,1]=PM_AQ_A0.integrate()
					
					PM_AQ_A1.setArgument(1, x[nn])
					A[1,2]=PM_AQ_A1.integrate()
					
					if (PM_AQ_A0.errorcode()==0 & PM_AQ_A1.errorcode()==0 & A[1]>0 & abs(A[2]/A[1])<abs(x[nn])) {
						PM[nn,.]=x[nn] - sign(A[2])*max((abs(A[2]/A[1]),0))
					}
				}
			}
		}
		
		// Pareto
		if (p_num==4) {
			PM_AQ_A0.setEvaluator(&PM_AQ_A0_4_f())
			PM_AQ_A0.setArgument(2, p_par)
			PM_AQ_A0.setArgument(3, cutoff)
		
			PM_AQ_A1.setEvaluator(&PM_AQ_A1_4_f())
			PM_AQ_A1.setArgument(2, p_par)
			PM_AQ_A1.setArgument(3, cutoff)
			
			PM_AQ_A0.setLimits((0, .))
			PM_AQ_A1.setLimits((0, .))
			for (nn=1; nn<=dim_x; nn++) {
				if (x[nn]==0) PM[nn,.]=0
				else {
					PM_AQ_A0.setArgument(1, x[nn])
					A[1,1]=PM_AQ_A0.integrate()
					
					PM_AQ_A1.setArgument(1, x[nn])
					A[1,2]=PM_AQ_A1.integrate()
					
					if (PM_AQ_A0.errorcode()==0 & PM_AQ_A1.errorcode()==0 & A[1]>0 & abs(A[2]/A[1])<abs(x[nn])) {
						PM[nn,.]=x[nn] - sign(A[2])*max((abs(A[2]/A[1]),0))
					}
				}
			}
		}
		
		// Cauchy
		if (p_num==5) {
			PM_AQ_A0.setEvaluator(&PM_AQ_A0_5_f())
			PM_AQ_A0.setArgument(2, cutoff)

			PM_AQ_A1.setEvaluator(&PM_AQ_A1_5_f())
			PM_AQ_A1.setArgument(2, cutoff)
			
			PM_AQ_A0.setLimits((0, .))
			PM_AQ_A1.setLimits((0, .))
			for (nn=1; nn<=dim_x; nn++) {
				if (x[nn]==0) PM[nn,.]=0
				else {
					PM_AQ_A0.setArgument(1, x[nn])
					A[1,1]=PM_AQ_A0.integrate()
					
					PM_AQ_A1.setArgument(1, x[nn])
					A[1,2]=PM_AQ_A1.integrate()
					
					if (PM_AQ_A0.errorcode()==0 & PM_AQ_A1.errorcode()==0 & A[1]>0 & abs(A[2]/A[1])<abs(x[nn])) {
						PM[nn,.]=x[nn] - sign(A[2])*max((abs(A[2]/A[1]),0))
					}
				}
			}
		}
		
		// Log
		if (p_num==6) {
			PM_AQ_A0.setEvaluator(&PM_AQ_A0_6_f())
			PM_AQ_A0.setArgument(2, p_par)
			PM_AQ_A0.setArgument(3, cutoff)

			PM_AQ_A1.setEvaluator(&PM_AQ_A1_6_f())
			PM_AQ_A1.setArgument(2, p_par)
			PM_AQ_A1.setArgument(3, cutoff)
			
			PM_AQ_A0.setLimits((0, .))
			PM_AQ_A1.setLimits((0, .))
			for (nn=1; nn<=dim_x; nn++) {
				if (x[nn]==0) PM[nn,.]=0
				else {
					PM_AQ_A0.setArgument(1, x[nn])
					A[1,1]=PM_AQ_A0.integrate()
					
					PM_AQ_A1.setArgument(1, x[nn])
					A[1,2]=PM_AQ_A1.integrate()
					
					if (PM_AQ_A0.errorcode()==0 & PM_AQ_A1.errorcode()==0 & A[1]>0 & abs(A[2]/A[1])<abs(x[nn])) {
						PM[nn,.]=x[nn] - sign(A[2])*max((abs(A[2]/A[1]),0))
					}
				}
			}
		}

		// Horseshoe
		if (p_num==7) {
			PM_AQ_A0.setEvaluator(&PM_AQ_A0_7_f())
			PM_AQ_A0.setArgument(2, p_par)
			PM_AQ_A0.setArgument(3, cutoff)

			PM_AQ_A1.setEvaluator(&PM_AQ_A1_7_f())
			PM_AQ_A1.setArgument(2, p_par)
			PM_AQ_A1.setArgument(3, cutoff)
			
			PM_AQ_A0.setLimits((0, 1))
			PM_AQ_A1.setLimits((0, 1))
			for (nn=1; nn<=dim_x; nn++) {
				if (x[nn]==0|abs(x[nn])>=600) PM[nn,.]=x[nn]
				else{
					PM_AQ_A0.setArgument(1, x[nn])
					A[1,1]=PM_AQ_A0.integrate()
					
					PM_AQ_A1.setArgument(1, x[nn])
					A[1,2]=PM_AQ_A1.integrate()			
					
					if (PM_AQ_A0.errorcode()==0 & PM_AQ_A1.errorcode()==0 & A[2]/A[1]<1 & A[2]/A[1]>0) {
						PM[nn,.]=(1-A[2]/A[1])*x[nn]
					}			
				}
			}
		}
	}
	
	// Return results on the posterior mean
	return(PM) 
}
//-----------------------------------------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------------------------------	
// Routine			: PLUG_IN
// Description		: Plug-in estimators of the bias and variance of the posterior mean 
//-----------------------------------------------------------------------------------------------------------------------	
real matrix PLUG_IN(real colvector hat_theta, real matrix tabmom, real scalar m_num, 		 		/*
		*/ 		real scalar p_num, real colvector p_par)
{
	// Determine the estimation method for each component of hat_theta
	step=0.01
	MT=colmax(tabmom[.,1])
	method=(abs(hat_theta):<=MT-step) +2:*(abs(hat_theta):>MT-step)
	T1=select(hat_theta, method:==1)
	T2=select(hat_theta, method:==2)
	m1_ind=selectindex(method:==1)
	m2_ind=selectindex(method:==2)

	// Method 1: linear interpolation between tabulaled moments (components of hat_theta \le 30)
	if (rows(T1)>0) {
		C=1/step
		T1_sign=sign(T1)
	
		T1_L =abs(trunc(T1:*C):/C)
		T1_U =abs(trunc((T1:*C):+T1_sign):/C)
		
		T1_L_ind=round(T1_L:*C:+1)
		T1_U_ind=round(T1_U:*C:+1)
		
		T1_L_wgt =editmissing((T1_U :-abs(T1)):/(T1_U :-T1_L),0)
		T1_U_wgt =(1:-T1_L_wgt)

		if (m_num==1) {
			hat_delta_1	=T1_sign:*(T1_L_wgt:*tabmom[T1_L_ind,2]+T1_U_wgt:*tabmom[T1_U_ind,2])
		}
		if (m_num==2) {
			hat_sigma2_1=T1_L_wgt :* tabmom[T1_L_ind,3] :+ T1_U_wgt:*tabmom[T1_U_ind,3]
		}
	}

	// Method 2: asymptotic approximations (components of hat_theta > 30)
	if (rows(T2)>0) {
		T2_sign=sign(T2)
		tab_points=rows(tabmom[.,1]) 
		//MC_theta_last=tabmom[tab_points,1]
		
		if (p_num==1) {
			pp_c		= p_par[1]
			psi			= J(rows(T2),1, pp_c)	
			bar_psi		= pp_c
			psi_d1		= J(rows(T2),1,0)	
			bar_psi_d1	= 0
		}
		if (p_num==2) { 
			pp_b		= p_par[1]
			pp_c		= p_par[2]
			psi			= pp_b:*pp_c:*abs(T2):^(pp_b-1) 
			bar_psi		= pp_b *pp_c *      MT^(pp_b-1)
			psi_d1		= pp_b:*pp_c:*(pp_b-1):*abs(T2):^(pp_b-2) 
			bar_psi_d1	= pp_b *pp_c *(pp_b-1) *      MT^(pp_b-2)
		}
		if (p_num==3) { 
			pp_b		= p_par[1]
			pp_c		= p_par[2]
			psi			= (1:-pp_b):/abs(T2):+pp_b:*pp_c:*abs(T2):^(pp_b-1) 					
			bar_psi		= (1 -pp_b) /MT      +pp_b *pp_c *      MT^(pp_b-1)
			psi_d1		=-(1:-pp_b):/(T2:^2):-pp_b:*pp_c:*(1:-pp_b):*abs(T2):^(pp_b-2) 		
			bar_psi_d1	=-(1 -pp_b) /( MT^2) -pp_b *pp_c *(1 -pp_b) *      MT^(pp_b-2)
		}
		if (p_num==4) { 
			pp_a		= p_par[1]
			pp_c		= p_par[2]
			psi			= pp_c:/(pp_a:*(1:+pp_c:*abs(T2)))
			bar_psi		= pp_c /(pp_a *(1 +pp_c *     MT))	
			psi_d1		=-pp_c^2:/(pp_a:*(1:+pp_c:*abs(T2)):^2)
			bar_psi_d1	=-pp_c^2 /(pp_a *(1 +pp_c *     MT)^2 )
		}
		if (p_num==5) { 
			psi			=2:*abs(T2):/(1:+(T2):^2)				
			bar_psi		=2:*     MT /(1 +   MT^2)				
			psi_d1		=2:*(1:-(T2):^2):/(1:+(T2):^2):^2		
			bar_psi_d1	=2 *(1 -   MT^2) /(1 +   MT^2)^2		
		}
		if (p_num==6) { 
			pp_c		=p_par[1]
			psi			=2:*pp_c:^2:/((abs(T2):^3:+pp_c:^2:*abs(T2)):*log(1:+pp_c:^2:/(T2:^2)))
			bar_psi		=2 * pp_c^2 /((      MT^3 +pp_c^2  *     MT) *log(1 +pp_c^2  /( MT^2)))
			psi_d1		=(4:*pp_c:^4:-2:*pp_c:^2:*(3:*T2:^2:+pp_c:^2):*log(1:+pp_c:^2:/(T2:^2))):/(((abs(T2):^3:+pp_c:^2:*abs(T2)):*log(1:+pp_c:^2:/(T2:^2))):^2)	
			bar_psi_d1	=(4*  pp_c^4 -2 * pp_c^2 *(3 * MT^2 + pp_c^2) *log(1 + pp_c^2 /( MT^2))) /(((      MT^3 + pp_c^2 *     MT) *log(1 + pp_c^2 /( MT^2)))^2 )	
		}
		if (p_num==7) {
			pp_c		= p_par[1]
			s			= (  T2:^2):/(2*pp_c^2)	
			bar_s		= (MT^2) /(2*pp_c^2)
			I_s 		= I_INT_AQ1(s)
			bar_I_s		= I_INT_AQ1(bar_s)
			psi			= 2:/(abs(T2):*     I_s):-(abs(T2):/(pp_c^2)) 
			bar_psi		= 2 /(     MT * bar_I_s) -(     MT /(pp_c^2))
			psi_d1		= (4:-4 :*     s:*     I_s :-2:*     I_s):/( (T2:^2) :*    (I_s:^2)):-pp_c^(-2) 
			bar_psi_d1	= (4 -4  * bar_s * bar_I_s  -2 * bar_I_s) /(  (MT^2)  * (bar_I_s^2)) -pp_c^(-2) 
		}			
		if (m_num==1) {
			bar_b=tabmom[tab_points,2]
			constant_bias	=MT^2*(1+bar_b/bar_psi)
			hat_delta_2		=-sign(T2):*(psi:*(1:-constant_bias:/(T2:^2))) 	
		}			
		if (m_num==2) {
			bar_v=tabmom[tab_points,3]
			constant_var	=(bar_v-1+2*bar_psi_d1)*(MT^2)/bar_psi
			hat_sigma2_2	=1:-2:*psi_d1:+constant_var:*psi:/(T2:^2)		
		}			
	}

	// Combine results from the two methods and return estimates
	if (m_num==1) {
		if (rows(T1)>0  & rows(T2)==0) hat_delta=hat_delta_1
		if (rows(T1)==0 & rows(T2)>0 ) hat_delta=hat_delta_2
		if (rows(T1)>0  & rows(T2)>0 ) hat_delta=sort((m1_ind,hat_delta_1)\(m2_ind,hat_delta_2),1)[.,2]
		return(hat_delta)
	}
	if (m_num==2) {
		if (rows(T1)>0  & rows(T2)==0) hat_sigma2=hat_sigma2_1
		if (rows(T1)==0 & rows(T2)>0 ) hat_sigma2=hat_sigma2_2
		if (rows(T1)>0  & rows(T2)>0 ) hat_sigma2=sort((m1_ind,hat_sigma2_1)\(m2_ind,hat_sigma2_2),1)[.,2]
		return(hat_sigma2)
	}
}
//-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// Routine			: DRAWS_BCPM_GL
// Description		: Draws of the (bias-corrected) posterior mean computed by Gauss-Laguerre/Gauss-Legendre
//-----------------------------------------------------------------------------------------------------------------------	
real matrix DRAWS_BCPM_GL(real colvector hat_theta, real scalar reps, 									/*
		*/ 		real scalar p_num, real colvector p_par, 												/*
		*/		real scalar NQP, real colvector QP, real colvector QW, real scalar cutoff,				/*
		*/		real scalar COD_BIAS, real matrix tabmom)
{
	K2=rows(hat_theta)
	draws=J(reps,K2,0)
	if (COD_BIAS==1) {
		for (h=1; h<=K2; h++) {
			xs_h=rnormal(reps, 1, hat_theta[h], 1)
			draws[.,h]=NLM_PM_GL(xs_h,p_num,p_par,NQP,QP,QW,cutoff) 
			draws[.,h]=draws[.,h]-PLUG_IN(draws[.,h],tabmom,1,p_num,p_par)     
		}
	}
	if (COD_BIAS==2) {
		for (h=1; h<=K2; h++) {
			xs_h=rnormal(reps, 1, hat_theta[h], 1)
			draws[.,h]=NLM_PM_GL(xs_h,p_num,p_par,NQP,QP,QW,cutoff) 
			draws[.,h]=draws[.,h]-PLUG_IN(xs_h,tabmom,1,p_num,p_par)     
		}		
	}		
	return(draws)
}
//-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// Routine			: DRAWS_BCPM_AQ
// Description		: Draws of the (bias-corrected) posterior mean computed by adaptive quadrature
//-----------------------------------------------------------------------------------------------------------------------	
real matrix DRAWS_BCPM_AQ(real colvector hat_theta, real scalar reps, 				/*
		*/ 		real scalar p_num, real colvector p_par, 							/*
		*/		real scalar AT, real scalar RT, real scalar cutoff,					/*
		*/		real scalar COD_BIAS, real matrix tabmom)
{
	K2=rows(hat_theta)
	draws=J(reps,K2,0)
	if (COD_BIAS==1) {
		for (h=1; h<=K2; h++) {
			xs_h=rnormal(reps,1,hat_theta[h],1)
			draws[.,h]=NLM_PM_AQ(xs_h,p_num,p_par,AT,RT,cutoff)	
			draws[.,h]=draws[.,h]-PLUG_IN(draws[.,h],tabmom,1,p_num,p_par) 
		}
	}
	if (COD_BIAS==2) {
		for (h=1; h<=K2; h++) {
			xs_h=rnormal(reps,1,hat_theta[h],1)
			draws[.,h]=NLM_PM_AQ(xs_h,p_num,p_par,AT,RT,cutoff)	
			draws[.,h]=draws[.,h]-PLUG_IN(xs_h,tabmom,1,p_num,p_par)	
		}
	}		
	return(draws)
}
//-----------------------------------------------------------------------------------------------------------------------	



end
//-----------------------------------------------------------------------------------------------------------------------	


//-----------------------------------------------------------------------------------------------------------------------	
// Di_wals: display of wals output 
//-----------------------------------------------------------------------------------------------------------------------	
program Di_wals
	syntax [, LEVel(cilevel) noTABle noFOCus noAUXiliary noHEADer cformat(string)]

	local y=abbrev(`"`e(depvar)'"',10)
	local focvars 	"`e(focvars)'"
	local auxvars	"`e(auxvars)'"
	local k1=e(k1)
	local k2=e(k2)
	local regressors: coln 	 e(b)
	local LSU_k		: colsof e(b) 
	local cons "`e(constype)'"
	if "`cons'"=="focus" 	{
		local k1_nc=`k1'-1
		local k2_nc=`k2'
		local focvars_nc: subinstr local focvars "_cons" ""
		local M_names "`focvars_nc' `auxvars' _cons"
	}
	else if "`cons'"=="auxiliary" {
		local k1_nc=`k1'
		local k2_nc=`k2'-1
		local M_names "`focvars' `auxvars'"
	}
	else {
		local k1_nc=`k1'
		local k2_nc=`k2'
		local M_names "`focvars' `auxvars'"
	}	
	tempname Bias RMSE 
	matrix `Bias'=e(bias)
	matrix `RMSE'=e(rmse)
	
	if "`header'"=="noheader" 		local hdis 0
	else 							local hdis 1
	if "`focus'"=="nofocus" 		local fdis 0
	else							local fdis 1
	if "`auxiliary'"=="noauxiliary" local adis 0
	else							local adis 1
	if "`table'"=="notable" 		local tdis 0
	else							local tdis 1

	* Check display format 
	if "`cformat'"=="" {
		local cformat "%8.0g"
	}
	else {
		local cfor_len: subinstr local cformat "%" " "
		local cfor_len: subinstr local cfor_len "." " "
		local cfor_len: word 1 of `cfor_len'
		cap confirm integer number `cfor_len'
		if _rc!=0 {
			local cformat "%8.0g"
		}
		else {
			if abs(`cfor_len')>8 {
				local cformat "%8.0g"
			}
		}	
	}
	
	local noncoll_X ""
	foreach jj of local regressors {
		_ms_parse_parts `jj'
		if `r(omit)'==0	local noncoll_X "`noncoll_X' `jj'"
	}
	local LL=10
	foreach xx of local noncoll_X {
	    local ll: strlen local xx
		local LL=max(`ll', `LL')
	}

	tempname CI WALS_CI 
	walsci `"`e(bcsimdata)'"' `level' `CI'
	matrix `WALS_CI'=J(2,`LSU_k',0)
	local ss=1
	forvalues jj=1(1)`LSU_k'{
		local xx: word `jj' of `regressors'
		local yy: word `ss' of `M_names'
		if "`xx'"=="`yy'" {
			matrix `WALS_CI'[1,`jj']=`CI'[1,`ss']
			matrix `WALS_CI'[2,`jj']=`CI'[2,`ss']
			local ++ss
		}
	}
	mata: st_replacematrix("e(ci)",st_matrix("`WALS_CI'"))

	local COEF_tit="      "
	if "`e(plugin)'"=="ds" {
		local MOM_tit	"  DS   "
		local CI_tit 	"  MC-DS   "
	}
	else {
		local MOM_tit	"  ML   "
		local CI_tit 	"  MC-ML   "
	}

	* Header
	if `hdis'==1 {
		local p1 =53
		local p2 =67
		#delimit;
		di as text _n `"`e(title)'"'					
			   as text _col(`p1') "Number of obs" 	_col(`p2') "="
			   as res %8.0f e(N)
			
			_n as text "Prior : " as res "`e(prior)'"
			   as text 	_col(`p1') "Residual df" 	_col(`p2') "=" 	
			   as res %8.0f e(df_r) 							
			
			_n as text _col(`p1') "k1 " 			_col(`p2') "=" 	
			   as res %8.0f e(k1) 							
			
			_n as text _col(`p1') "k2 " 			_col(`p2') "=" 	
			   as res %8.0f e(k2) 									
			
			_n as text _col(`p1') "sigma " 			_col(`p2') "=" 
			   as res %8.0g e(sigma) 								
			
			_n as text _col(`p1') "MC reps " 		_col(`p2') "=" 
			   as res %8.0f e(reps); 								
		#delimit cr	
	}
	
	* Coefficient table
	if `tdis'==1 {
		local LLp1=`LL'+1
		local p2=`LL'+8	
		local p3=`LL'+18
		local p4=`LL'+28
		local p5=`LL'+38
		local p6=`LL'+53
		#delimit;
		di as text "{hline `LLp1'}{c TT}{hline 63}"													
			_n _col(`LL') "  {c |}"																
			_col(`p2') "`COEF_tit'"                            									
			_col(`p3') "`MOM_tit'"                            									
			_col(`p4') "`MOM_tit'"														 
			_col(`p5') "`MOM_tit'"														 
			_col(`p6') "`CI_tit'" 	                           									
			_n "{ralign `LL':`y'} {c |}"															
			"     Coef.      Bias    Std.Err.    RMSE      [`level'% Conf. Int.]" 
			_n "{hline `LLp1'}{c +}{hline 63}";	
		#delimit cr	

		local p1=`LL'+3	
		local p2=`LL'+16	
		local p3=`LL'+26
		local p4=`LL'+35
		local p5=`LL'+47
		local p6=`LL'+58
		local hh=1	
		local xx_num=1
		foreach xx of local regressors {
			_ms_parse_parts `xx'
			if `r(omit)'==0 {
				if (`hh'==1 & `k1_nc'>0 & `fdis'==1) {
				    di as text "focus" _col(`LL') "  {c |}"
				}
				if (`hh'==`k1_nc'+1 & `k1_nc'>0 & `k2_nc'>0 & `adis'==1)	{
					if `fdis'==1	di as text "{hline `LLp1'}{c +}{hline 63}"
					di as text "auxiliary" _col(`LL') "  {c |}"
				}
				if (`hh'==1 & `k1_nc'==0 & `k2_nc'>0 & `adis'==1) {
					di as text "auxiliary" _col(`LL') "  {c |}"
				}
				if `hh'==`k1_nc'+`k2_nc'+1 & "`cons'"!="noconstant" {
					if "`cons'"=="focus" & `fdis'==1 {
					    if `k2_nc'>0 di as text "{hline `LLp1'}{c +}{hline 63}"
						di as text "focus" _col(`LL') "  {c |}"
					}
					if "`cons'"=="auxiliary" & `adis'==1 {
					    if `k2_nc'>0 di as text "{hline `LLp1'}{c +}{hline 63}"
						di as text "auxiliary" _col(`LL') "  {c |}"
					}
				}
				if 		(`hh'<=`k1_nc' & `fdis'==1)											///
					|	(`hh'>=`k1_nc'+1 & `hh'<=`k1_nc'+`k2_nc' & `k2_nc'>0 & `adis'==1)	///
					|	(`hh'==`k1_nc'+`k2_nc'+1 & "`cons'"=="focus" & `fdis'==1)			///
					|	(`hh'==`k1_nc'+`k2_nc'+1 & "`cons'"=="auxiliary" & `adis'==1)		{
					
					#delimit;
					di as text "{ralign `LL':`xx'} {c |}  "						
							_col(`p1') as res `cformat'  _b[`xx'] 				"  "
							_col(`p2') as res `cformat'  `Bias'[1,`xx_num']		"  "
							_col(`p3') as res `cformat'  _se[`xx']				"  "
							_col(`p4') as res `cformat'  `RMSE'[1,`xx_num']		"  "
							_col(`p5') as res `cformat' `CI'[1,`hh'] 			"  "
							_col(`p6') as res `cformat' `CI'[2,`hh'];		
					#delimit cr	
				}
				local hh=`hh'+1
			}
			local xx_num=`xx_num'+1
		}
		noi di as text "{hline `LLp1'}{c BT}{hline 63}"		
	}
end
//-----------------------------------------------------------------------------------------------------------------------	


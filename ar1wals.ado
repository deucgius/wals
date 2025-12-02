* Weighted Average Least Squares (WALS) with AR(1) errors
*! Version 1.0
* Date			: 2024/10/05
* Authors		: De Luca Giuseppe, Jan R. Magnus
*-----------------------------------------------------------------------------------------------------------------------	
* ar1wals
*-----------------------------------------------------------------------------------------------------------------------	
program define ar1wals, eclass sort
	version 14.0, missing
	if replay() {
		if "`e(cmd)'" != "ar1wals" {
			error 301
		}
		syntax [, SAVing(string asis) noHEADer noTABle noFOCus noAUXiliary LEVel(cilevel) cformat(string)]
		local bcsimdata `"`e(bcsimdata)'"'
		cap confirm file `"`bcsimdata'"'
		if _rc {
			di as err "MC replicates of bias-corrected ar1wals estimates not found. You need to refit the model and resave the results" 
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
		Di_ar1wals, `header' `table' `focus' `auxiliary' lev(`level') `cformat'
		exit
	}
	if !replay() {
		syntax [anything] [if] [in], [*]
		local cmdline : copy local 0
	}
	Estimate `0'
	ereturn local cmdline `"ar1wals `0'"'
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
* Syntax (weights not allowed)
*------------------------------------------------------------------------------------------------
#delimit;	
syntax [anything] [if] [in], 											
	[														
	AUXCONstant												
	noCONstant												
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
	
	CORC
	RHOtype(string)
	SSEsearch
	TWOstep
	ITERate(integer 300) 
	TOLerance(real 1e-6)
	noLOg
	showar1
	
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
	exit 198
}
if `:list depvar in auxlist0' {
	di as err "dependent variable {bf:`depvar'} not allowed as auxiliary regressor"
	exit 198
}
if subinstr("`foclist0' `auxlist0'",".`depvar' ","",.) != "`foclist0' `auxlist0'" {
	di as err "time-series operators of dependent variable {bf:`depvar'} may not be " /*
	*/ "included as independent variables"
	exit 198
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Mark estimation sample 
*------------------------------------------------------------------------------------------------
marksample touse
_ts tvar
markout `touse' `depvar' `foclist0' `auxlist0' `tvar'
count if `touse' 
if r(N)==0 	error 2000 
if r(N)==1 	error 2001
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Preliminary check of options 
*------------------------------------------------------------------------------------------------
if "`rhotype'"=="" 	local rhotype 	 "regress"
if "`showar1'"=="" 	local show_prais "qui"
else 				local show_prais "noi"
local prais_opts "`corc' `twostep' `ssesearch' rhotype(`rhotype') iterate(`iterate') tolerance(`tolerance') `log'"

* saving
cap erase `"`c(tmpdir)'/STWALS_000001.tmp"'
noi _savingopt_parse CIMC_save CIMC_savereplace : saving ".dta" `"`saving'"'
if `"`CIMC_save'"'!="" & `"`CIMC_savereplace'"'=="" confirm new file `"`CIMC_save'"'
if `"`CIMC_save'"'=="" {
	local CIMC_save `"`c(tmpdir)'/STWALS_000001.tmp"'
	local CIMC_savereplace "replace"
}

* display options
if "`fast'"=="" {
	if "`cformat'"!="" local cformat "cformat(`cformat')"
	local display_options "`table' `focus' `auxiliary' `header' `cformat'"
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* prior (default pareto)
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
* quadmethod and quadext
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
}
if "`quadm'"=="adaptive" {		
	if "`quadext'"!=""		noi di as text "note: {bf:quadext} is ineffective with {bf:quadmethod(adaptive)}"
	if "`quadnpts'"!="500"	noi di as text "note: {bf:quadnpts} is ineffective with {bf:quadmethod(adaptive)}"
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Constant term (if any)
*------------------------------------------------------------------------------------------------
local nc `constant'
local ac "`auxconstant'"
if "`nc'"~="noconstant" {
	tempvar const
	gen double `const'=1 if `touse'
}
else {
	if "`ac'"!="" {
		noi di in red "cannot specify both {bf:noconstant} and {bf:auxconstant}"
		error 184
	}
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* check collinearity and estimate rho
*------------------------------------------------------------------------------------------------
* Only for rmcoll message 
if "`showar1'"=="" {
	noi version 11: _rmdcoll `depvar' `foclist0' `auxlist0' if `touse', `nc'
}

* Prais estimates of unrestricted model
cap `show_prais' prais `depvar' `foclist0' `auxlist0' if `touse', `nc' `prais_opts'
if _rc!=0 {
	noi "Convergence of prais estimates for unrestricted model not achieved"
	error 430
}	
tempvar touse2
gen `touse2'=e(sample)
tempname rho dw dw_0
scalar `rho'=e(rho)
scalar `dw'=e(dw)
scalar `dw_0'=e(dw_0)
local LSU_k: colsof e(b)
local LSU_names: coln e(b)
local N=e(N)
local rho_method "`e(method)'" 
local rho_type	 "`e(rhotype)'"

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
local TEMP_foc_vars ""
foreach xx of local f_names {
	fvrevar `xx'
	local TEMP_foc_vars 	"`TEMP_foc_vars' `r(varlist)'"
}
if "`nc'"~="noconstant" & "`ac'"=="" 	local TEMP_foc_vars `const' `TEMP_foc_vars' 
local TEMP_foc_names 	`f_names'
if "`nc'"~="noconstant" & "`ac'"==""	local TEMP_foc_names _cons `TEMP_foc_names'

* List of noncollinear auxiliary regressors (as temp vars), including the constant term (if any)
local TEMP_aux_vars ""
foreach xx of local a_names {
	fvrevar `xx'
	local TEMP_aux_vars 	"`TEMP_aux_vars' `r(varlist)'"
}
if "`nc'"~="noconstant" & "`ac'"!="" 	local TEMP_aux_vars  `TEMP_aux_vars' `const' 
local TEMP_aux_names 	`a_names'
if "`nc'"~="noconstant" & "`ac'"!=""	local TEMP_aux_names `TEMP_aux_names' _cons

* List of all noncollinear regressors (as temp vars), including the constant term (if any)
local TEMP_names `TEMP_foc_names' `TEMP_aux_names'
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Number of focus and auxiliary regressors 
*------------------------------------------------------------------------------------------------
local k1: word count `TEMP_foc_vars' 
local k2: word count `TEMP_aux_vars' 
if `k2'<1 {
	noi di as err "The model must contain at least one auxiliary regressor (including, if any, the constant term)" 
	error 102
}
local k = `k1'+`k2'
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* AR(1) transform and WALS estimates 
*------------------------------------------------------------------------------------------------
* AR(1) transform
preserve 
keep if `touse'
local varlist `depvar' `TEMP_foc_vars' `TEMP_aux_vars' 
tempvar TVar
gen double `TVar'=`tvar'
tsset `TVar', noquery 
ar1_transf `varlist', rho(`rho') `corc' touse(`touse') 

* WALS estimates on transformed data
if "`TEMP_foc_vars'"!="" {
	#delimit;
	noi wals `depvar' (`TEMP_foc_vars') `TEMP_aux_vars'  
		if `touse2'==1, 
		nocons sigma(`sigma') prior(`PRI') `QUADRATURE' 
		plugin(`plugin') patht(`pathtab')
		lev(`level') rep(`reps') rseed(`rseed') saving(`saving') 
		nocollcheck notable nohead `fast';
	#delimit cr
}
else {
	#delimit;
	noi wals `depvar' `TEMP_aux_vars'  
		if `touse2'==1, 
		nocons sigma(`sigma') prior(`PRI') `QUADRATURE' 
		plugin(`plugin') patht(`pathtab')
		lev(`level') rep(`reps') rseed(`rseed') saving(`saving') 
		nocollcheck notable nohead `fast';
	#delimit cr
}	
if e(k1)!=`k1'|e(k2)!=`k2' {
	noi di as err "Collinear regressors in the transformed data"
	exit 198
}
if e(N)!=`N' {
	noi di as err "warning: differences in the sample size of prais and WALS estimates"
	exit 198
}
restore
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Parse wals results
*------------------------------------------------------------------------------------------------
tempname stddev prior_par M_W M_pm M_t_ratios
scalar `stddev'=e(sigma)
matrix `prior_par'=e(priorpar)
matrix `M_W'=e(wals_wgt)
matrix `M_pm'=e(wals_pm)
matrix `M_t_ratios'=e(t_ratios)
if "`fast'"=="" {
	local level=e(level)
	local reps=e(reps)
	local BIAS_type	"`e(plugin)'"
}
local df_r=e(df_r)
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Estimated coefficients and their moments 
*------------------------------------------------------------------------------------------------
tempname M0_b 
matrix `M0_b'		=e(b)
if "`fast'"=="" {
	tempname M0_bias M0_V M0_MSE M0_rmse  
	matrix `M0_bias'	=e(bias)
	matrix `M0_V'		=e(V)
	matrix `M0_MSE'		=e(MSE)
	matrix `M0_rmse'	=e(rmse)
}

* If any, put constant term as last element  
tempname M_b 
if "`nc'"~="noconstant" & "`ac'"=="" {
	local M_names: subinstr local TEMP_names "_cons" ""
	local M_names `M_names' _cons
	local M_foc_names: subinstr local TEMP_foc_names "_cons" ""
	local M_foc_names `M_foc_names' _cons
	local M_aux_names `TEMP_aux_names'
	matrix `M_b'	=(`M0_b'[1,2..`k'],`M0_b'[1,1])
	matrix coln `M_b' 	= `M_names'
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
    local M_names `TEMP_names'
	local M_foc_names `TEMP_foc_names'
	local M_aux_names `TEMP_aux_names'
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

* List of omitted regressors
local M_focaux_names `M_foc_names' `M_aux_names' 
local OMITTED_vars: list LSU_names - M_focaux_names
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Finalize dataset containing the MC draws of bias-corrected WALS estimates
*------------------------------------------------------------------------------------------------
if "`fast'"=="" {
	local CIMC_save `"`e(bcsimdata)'"'
	preserve
	use `"`CIMC_save'"', clear
	local jj=1
	foreach xx in `TEMP_foc_names' `TEMP_aux_names' {
		if "`nc'"~="noconstant" & "`ac'"=="" {
			if `jj'==1 	rename wals_bc_`jj' wals_bc_0
			else {
				local ll=`jj'-1
				rename wals_bc_`jj' wals_bc_`ll'
				lab var wals_bc_`ll' "`xx'" 
			}
		}
		else {
			rename wals_bc_`jj' wals_bc_`jj'
			lab var wals_bc_`jj' "`xx'" 
		}
		local jj=`jj'+1
	}
	if "`nc'"~="noconstant" & "`ac'"=="" {
		rename wals_bc_0 wals_bc_`k'
		lab var wals_bc_`k' "`xx'" 
		order wals_bc_`k', last
	}
	saveold `"`CIMC_save'"', `CIMC_savereplace'
	local wals_bc_fname `"`c(filename)'"'
	restore
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Return estimation results 
*------------------------------------------------------------------------------------------------
if "`fast'"!="" {
	ereturn post `WALS_b' , dep(`depvar') obs(`N') esample(`touse2') dof(`df_r') buildfvinfo 
	ereturn scalar k0 		=`LSU_k'
	ereturn scalar k1 		=`k1'
	ereturn scalar k2 		=`k2'
	ereturn scalar df_r		=`e(df_r)'				
	ereturn scalar sigma	=`stddev'
	ereturn scalar rho		=`rho'	
	ereturn scalar dw_0		=`dw_0'	
 	ereturn scalar dw		=`dw'	
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
	
	ereturn local predict 		"ar1wals_p"
	ereturn local properties	"`e(properties)'"
	if "`PRI'"=="laplace" 		ereturn local quadmethod	"analytic"
	else 						ereturn local quadmethod	"`quadm'"
	ereturn local prior 		"`PRI'"
	if "`corc'" == "" 	ereturn local tranmeth "prais"
	else 				ereturn local tranmeth "corc"
	ereturn local rhomethod 	"`rho_method'" 
	ereturn local rhotype		"`rho_type'"	
	if "`nc'"=="noconstant"						ereturn local constype "`nc'"
	else if "`nc'"~="noconstant" & "`ac'"=="" 	ereturn local constype "focus"
	else 										ereturn local constype "auxiliary"
	ereturn local auxvars 		"`M_aux_names'"
	ereturn local focvars 		"`M_foc_names'"
	ereturn local omitvars	 	"`OMITTED_vars'"
	ereturn local allvars 		"`LSU_names'"
	ereturn local depvar		"`e(depvar)'"
	ereturn local title 		"AR(1) WALS estimates"
	ereturn local cmd 			"ar1wals"
	
	ereturn matrix priorpar		=`prior_par'
 	ereturn matrix wals_wgt		=`M_W'
 	ereturn matrix wals_pm		=`M_pm'
 	ereturn matrix t_ratios		=`M_t_ratios'	
}
else{
	ereturn post `WALS_b' `WALS_V', dep(`depvar') obs(`N') esample(`touse2') dof(`df_r') buildfvinfo
 	ereturn scalar k0 		=`LSU_k'
 	ereturn scalar k1 		=`k1'
 	ereturn scalar k2 		=`k2'
 	ereturn scalar rank 	=`k'
	ereturn scalar df_r		=`e(df_r)'				
 	ereturn scalar sigma	=`stddev'
 	ereturn scalar rho		=`rho'	
	ereturn scalar dw_0		=`dw_0'	
 	ereturn scalar dw		=`dw'	
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
 	ereturn scalar reps 	=`reps'
 	ereturn scalar level    =`level'

	ereturn local marginsok 	"XB default"
	ereturn local predict 		"ar1wals_p"
	ereturn local properties	"`e(properties)'"
	ereturn local bcsimdata 	`"`wals_bc_fname'"'	
	ereturn local plugin		"`plugin'"
	if "`PRI'"=="laplace" 		ereturn local quadmethod	"analytic"
	else 						ereturn local quadmethod	"`quadm'"
	ereturn local prior 		"`PRI'"
	if "`corc'" == "" 	ereturn local tranmeth "prais"
	else 				ereturn local tranmeth "corc"
	ereturn local rhomethod 	"`rho_method'" 
	ereturn local rhotype		"`rho_type'"	
	if "`nc'"=="noconstant"						ereturn local constype "`nc'"
	else if "`nc'"~="noconstant" & "`ac'"=="" 	ereturn local constype "focus"
	else 										ereturn local constype "auxiliary"
	ereturn local auxvars 		"`M_aux_names'"
	ereturn local focvars 		"`M_foc_names'"
	ereturn local omitvars	 	"`OMITTED_vars'"
	ereturn local allvars 		"`LSU_names'"
	ereturn local depvar		"`e(depvar)'"
	ereturn local title 		"AR(1) WALS estimates"
	ereturn local cmd 			"ar1wals"

	ereturn matrix priorpar		=`prior_par'
 	ereturn matrix wals_wgt		=`M_W'
 	ereturn matrix wals_pm		=`M_pm'
 	ereturn matrix t_ratios		=`M_t_ratios'
	ereturn matrix ci			=`WALS_CI'
	ereturn matrix MSE			=`WALS_MSE'
	ereturn matrix rmse			=`WALS_rmse'
	ereturn matrix bias			=`WALS_bias'
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Display estimation results 
*------------------------------------------------------------------------------------------------
if "`fast'"=="" 	noi Di_ar1wals, level(`level') `display_options'   
*------------------------------------------------------------------------------------------------
}
end
*-----------------------------------------------------------------------------------------------------------------------	


*-----------------------------------------------------------------------------------------------------------------------	
* AR(1) transform
*-----------------------------------------------------------------------------------------------------------------------	
program define ar1_transf
	syntax varlist, rho(string) touse(varname) [corc]
	tokenize `varlist'
	while ("`1'"!="") {
		tempvar tv
		_ms_parse_parts `1'
		if r(type) == "interaction" {
			local c = r(k_names)
			local l1
			local sharp
			forval j = 1/`c' {
				local l1 `l1'`sharp'`r(op`j')'l.`r(name`j')'
				local sharp "#"
			}
		}
		else {
			local l1 l.`1'
		}
		gen double `tv' = `1' - `rho'*`l1'  if `touse' & l.`touse'
		if "`corc'" == "" {
			replace `tv' = `1' * sqrt(1-`rho'^2)  if `touse' & l.`touse'!= 1
		}
		drop `1'
		rename `tv' `1'
		mac shift
	}
end	
*-----------------------------------------------------------------------------------------------------------------------	


//-----------------------------------------------------------------------------------------------------------------------	
// Di_ar1wals: display of wals output 
//-----------------------------------------------------------------------------------------------------------------------	
program Di_ar1wals
	syntax [, LEVel(cilevel) noTABle noFOCus noAUXiliary noHEADer cformat(string)]

	local y=abbrev(`"`e(depvar)'"',10)
	local focvars 	"`e(focvars)'"
	local auxvars	"`e(auxvars)'"
	local k1=e(k1)
	local k2=e(k2)
	local regressors: coln 	 e(b)
	local LSU_k		: colsof e(b) 
	if "`e(tranmeth)'"=="corc" 	local tranf "Cochrane-Orcutt"
	else 						local tranf "Prais-Winsten"
	 
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
	local chi : di %8.2f e(het_chi2)
	local chi = trim("`chi'")
	
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
		if "`e(rhomethod)'"!="SSE search" {
			#delimit;
			di as text _n `"`e(title)'"'					
				   as text _col(`p1') "Number of obs" 	_col(`p2') "="
				   as res %8.0f e(N)
				
				_n as text "Prior         : " as res "`e(prior)'"
				   as text 	_col(`p1') "Residual df" 	_col(`p2') "=" 	
				   as res %8.0f e(df_r) 							
				
				_n as text "Transformation: " as res "`tranf'" 
				   as text _col(`p1') "k1 " 			_col(`p2') "=" 	
				   as res %8.0f e(k1) 							
				
				_n as text "rho method    : " as res "`e(rhomethod)'"
				   as text _col(`p1') "k2 " 			_col(`p2') "=" 	
				   as res %8.0f e(k2) 									
				
				_n as text "rho type      : " as res "`e(rhotype)'"
				   as text _col(`p1') "sigma " 			_col(`p2') "=" 
				   as res %8.0g e(sigma) 								
				
				_n as text "rho = " as res %8.0g e(rho)
				   as text _col(`p1') "MC reps " 		_col(`p2') "=" 
				   as res %8.0f e(reps); 								
			#delimit cr	
		}
		else {
			#delimit;
			di as text _n `"`e(title)'"'					
				   as text _col(`p1') "Number of obs" 	_col(`p2') "="
				   as res %8.0f e(N)
				
				_n as text "Prior         : " as res "`e(prior)'"
				   as text 	_col(`p1') "Residual df" 	_col(`p2') "=" 	
				   as res %8.0f e(df_r) 							
				
				_n as text "Transformation: " as res "`tranf'" 
				   as text _col(`p1') "k1 " 			_col(`p2') "=" 	
				   as res %8.0f e(k1) 							
				
				_n as text "rho method    : " as res "`e(rhomethod)'"
				   as text _col(`p1') "k2 " 			_col(`p2') "=" 	
				   as res %8.0f e(k2) 									
				
				_n as text "rho = " as res %8.0g e(rho)
				   as text _col(`p1') "sigma " 			_col(`p2') "=" 
				   as res %8.0g e(sigma) 								
				
				_n as text _col(`p1') "MC reps " 		_col(`p2') "=" 
				   as res %8.0f e(reps); 								
			#delimit cr	
		}
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


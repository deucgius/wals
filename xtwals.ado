* Weighted Average Least Squares (WALS) for panel data
*! Version 1.0
* Date			: 2024/10/05
* Authors		: De Luca Giuseppe, Jan R. Magnus
*-------------------------------------------------------------------------------------------------
* xtwals
*-------------------------------------------------------------------------------------------------
program define xtwals, eclass sort
	version 14.0, missing
	if replay() {
		if "`e(cmd)'" != "xtwals" {
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
		Di_xtwals, `header' `table' `focus' `auxiliary' lev(`level') `cformat'
		exit
	}
	if !replay() {
		syntax [anything] [fw aw iw /] [if] [in], [*]
		local cmdline : copy local 0
	}
	Estimate `0'
	ereturn local cmdline `"xtwals `0'"'
end
*-------------------------------------------------------------------------------------------------



*-------------------------------------------------------------------------------------------------
* Estimate
*-------------------------------------------------------------------------------------------------
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
	/* noCONstant	- not allowed */
	/* sigma(real 0) - not allowed */ 											

	PRIor(string) 	
	
	QUADMethod(string)					
	QUADNPts(integer 500)		
	QUADExt(string)					
	QUADATol(real 1e-9) 
	QUADRTol(real 1e-7)		
	
	PLUGin(string)											
	PATHTab(string)
	
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
	
	RE
	ar1

	sa
	RHOType(string)
	RHOF(real -2)
	TWOstep
	SHOWXTreg
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
	*/ 	"included as independent variables"
	exit 198
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Verify the integrity of the panel data and mark estimation sample 
*------------------------------------------------------------------------------------------------
marksample touse
if "`ar1'"=="" {
	_xt
	local ivar "`r(ivar)'"
	markout `touse' `depvar' `foclist0' `auxlist0' `exp' `ivar'
}
else {
	_xt, trequired
	local ivar "`r(ivar)'"
	local tvar "`r(tvar)'"
	tempname tdelta
	scalar `tdelta' = r(tdelta)
	markout `touse' `depvar' `foclist0' `auxlist0' `exp' `ivar' `tvar' 
}	 			
count if `touse' 
if r(N)==0 	error 2000 
if r(N)==1 	error 2001
if subinword("`foclist0' `auxlist0'","`ivar'","",.) != "`foclist0' `auxlist0'" {
	di as err "panel variable {bf:`ivar'} may not be included as an independent variable" 
	exit 198
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Constant term 
*------------------------------------------------------------------------------------------------
local ac "`auxconstant'"
tempvar const
gen double `const'=1 if `touse'
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Preliminary check of options
*------------------------------------------------------------------------------------------------
if "`showxtreg'"!="" 	local show "noi"
if "`rhof'"!="-2" 		local rhof_ops "rhof(`rhof')"
if "`rhotype'"!=""		local rhot_ops "rhotype(`rhotype')"
if "`rhof_ops'"!="" & "`rhot_ops'"!="" {
		noi di as err "cannot specify both {bf:rhof} and {bf: rhotype}"
		error 198
}
if "`rhof_ops'"!="" & "`twostep'"!="" {
		noi di as err "cannot specify both {bf:rhof} and {bf: twostep}"
		error 198
}
if "`ar1'"=="" {
	if "`rhot_ops'"!="" {
		di as err "{bf:`rhot_ops'} not allowed without {bf:ar1}"
		error 198
	}
	if "`rhof_ops'"!="" {
		di as err "{bf:`rhof_ops'} not allowed without {bf:ar1}"
		error 198
	}
	if "`twostep'"!="" {
		di as err "{bf:twostep} not allowed without {bf:ar1}"
		error 198
	}
}
if "`re'"=="" {
	if "`ac'"!="" {
		noi di as err "{bf:auxconstant} not allowed with {bf:fe}"
		error 198
	}
	if "`sa'"!="" {
		di as err "{bf:sa} not allowed with {bf:fe}"
		error 198
	}
	if "`ar1'"=="" {
		local xtcomm "`show' xtreg"
		local xtcomm_opts "fe"
	}
	else {
		local xtcomm "`show' xtregar"
		local xtcomm_opts "fe `rhot_ops' `rhof_ops' `twostep'"
	}
}
if "`re'"!="" {
	if "`ar1'"=="" {
		local xtcomm "`show' xtreg"
		local xtcomm_opts "re `sa'"
	}
	else {
		if "`sa'"!="" {
			di as err "option {bf:sa} not allowed with {bf:ar1}"
			error 198
		}
		local xtcomm "`show' xtregar"
		local xtcomm_opts "re `rhot_ops' `rhof_ops' `twostep'"
	}
}

* Saving
cap erase `"`c(tmpdir)'/STWALS_000001.tmp"'
noi _savingopt_parse CIMC_save CIMC_savereplace : saving ".dta" `"`saving'"'
if `"`CIMC_save'"'!="" & `"`CIMC_savereplace'"'=="" confirm new file `"`CIMC_save'"'
if `"`CIMC_save'"'=="" {
	local CIMC_save `"`c(tmpdir)'/STWALS_000001.tmp"'
	local CIMC_savereplace "replace"
}

* Display options
if "`fast'"=="" {
	if "`cformat'"!="" local cformat "cformat(`cformat')"
	local display_options "`table' `focus' `auxiliary' `header' `cformat'"
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Weights 
*------------------------------------------------------------------------------------------------
if `"`weight'"' != "" {
	tempvar wgt_var	
	gen double `wgt_var'=`exp' 					if `touse'
	xtsum `wgt_var'								
	if abs(r(sd_w)) > 1e-8 {
		di as err "weights must be constant within panel"
		exit 198
	}
	if "`re'" != "" {
		noi di "weights not allowed with random-effects" 
		error 198
	}
	if (`"`weight'"'=="iweight"|`"`weight'"'=="pweight") {
		noi di as err "`weight' not allowed"
		error 198
	}
	else if `"`weight'"'=="aweight" {
		sum `wgt_var' , meanonly
		replace `wgt_var' = `wgt_var'/r(mean)	if `touse'
		local wgt_exp `"[aw=`wgt_var']"'
		local wgt_exp_post `"[aw=`exp']"'		
	}
	else {
		local wgt_exp `"[fw=`wgt_var']"'
		local wgt_exp_post `"[fw=`exp']"'
	}
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
			noi di in red "mata matrix {bf:quadext(`quadext')} not found"
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
* Check for perfect collinearity in the unrestricted model
*------------------------------------------------------------------------------------------------
* rmcoll message 
if "`showxtreg'"=="" {
	if "`re'"=="" {
		noi _regress `depvar' `foclist0' `auxlist0' `wgt_exp' if `touse', absorb(`ivar') nohead notable
	}
	else {
		noi _regress `depvar' `foclist0' `auxlist0' `wgt_exp' if `touse', nohead notable
	}
}

* LS estimates of unrestricted model	
`xtcomm' `depvar' `foclist0' `auxlist0' `wgt_exp' if `touse', `xtcomm_opts' 
local LSU_k: colsof e(b)
local LSU_names: coln e(b)

* Estimation sample from unrestricted model (only for the model without ar1)
if "`re'"!=""|"`ar1'"=="" {
	replace `touse'=e(sample)
}
if "`re'"!="" {
	tempname sig_e sig_u 
	scalar `sig_e'=e(sigma_e)
	scalar `sig_u'=e(sigma_u)
}
if "`ar1'"!="" {
	tempname rho
	scalar `rho'=e(rho_ar)
	local rho_type	 "`e(rhotype)'"
}

* List of noncollinear regressors, excluding the constant term
_ms_omit_info e(b)
tempname x0_omit
matrix `x0_omit'=r(omit)
local x_names ""
local jj=1
foreach xx of local LSU_names {
	*if `x0_omit'[1,`jj']==0 & "`xx'"!="_cons"	local x_names "`x_names' `xx'"
	if "`ar1'"=="" & `x0_omit'[1,`jj']==0 & "`xx'"!="_cons"						local x_names "`x_names' `xx'"
	if "`ar1'"!="" & `x0_omit'[1,`jj']==0 & "`xx'"!="_cons" & e(b)[1,`jj']!=0	local x_names "`x_names' `xx'"
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
if "`ac'"==""	local TEMP_foc_vars `const' `TEMP_foc_vars' 
local TEMP_foc_names 	`f_names'
if "`ac'"==""	local TEMP_foc_names _cons `TEMP_foc_names'

* List of noncollinear auxiliary regressors (as temp vars), including the constant term (if any)
local TEMP_aux_vars ""
foreach xx of local a_names {
	fvrevar `xx'
	local TEMP_aux_vars 	"`TEMP_aux_vars' `r(varlist)'"
}
if "`ac'"!="" 	local TEMP_aux_vars  `TEMP_aux_vars' `const' 
local TEMP_aux_names 	`a_names'
if "`ac'"!=""	local TEMP_aux_names `TEMP_aux_names' _cons

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
* WALS on trasformed data
*------------------------------------------------------------------------------------------------
preserve 
keep if `touse'
local N_g=e(N_g)
local g_min=e(g_min)
local g_avg=e(g_avg)
local g_max=e(g_max) 
local T_con=e(Tcon)

* Fixed-effects 
if "`re'"=="" {
	local sigma_W=e(sigma_e)
	local df_r=e(df_r)
	if "`ar1'"=="" {						/* iid errors 	*/

		* List of variables to transform 
		local varlist `depvar' `TEMP_foc_vars' `TEMP_aux_vars'

		* FE transform
		sort `ivar' 
		fe_transf `varlist' `wgt_exp', ivar(`ivar')
	}
	else {									/* AR(1) errors */

		* List of variables to transform 
		local varlist `depvar' `TEMP_foc_vars' `TEMP_aux_vars'
		
		* Time differences  
		sort `ivar' `tvar'
		tempvar difft
		by `ivar': gen double `difft'=(`tvar'[_n]-`tvar'[_n-1])/`tdelta' if _n>1
		
		* Ci_rho transform 
		sort `ivar' `tvar'
		Ci_rho_transf `varlist', ivar(`ivar') dif_t(`difft') rho(`rho')

		* Drop first obs by panel
		sort `ivar' `tvar'
		by `ivar': drop if _n==1
		
		* FE transform 
		sort `ivar' `tvar'
		fe_transf `varlist' `wgt_exp', ivar(`ivar')
	}
}

* Random-effects 
if "`re'"!="" {
	local varlist `depvar' `TEMP_foc_vars' `TEMP_aux_vars'
	local sigma_W=0	
	if "`ar1'"=="" {						/* iid errors 	*/

		* theta_i
		if `T_con'==1 {
			tempname theta_i
			scalar `theta_i' = e(theta)
		}
		else {
			tempvar T_i theta_i
			bys `ivar': gen `T_i'=_N 
			gen double `theta_i'=1-`sig_e'/sqrt(`T_i'*`sig_u'^2+`sig_e'^2) 
			sum `theta_i', d
			tempname theta_sum
			matrix `theta_sum'=(r(min),r(p5),r(p50),r(p95),r(max))
			matrix coln `theta_sum' = min p5 p50 p95 max
		}
		
		* re transform
		sort `ivar' 
		re_transf `varlist', ivar(`ivar') thetai(`theta_i')		
	}
	else {									/* AR(1) errors */

		* Time differences  
		sort `ivar' `tvar'
		tempvar difft
		by `ivar': gen double `difft'=(`tvar'[_n]-`tvar'[_n-1])/`tdelta' if _n>1

		* Ci_rho transform 
		sort `ivar' `tvar'
		Ci_rho_transf `varlist', ivar(`ivar') dif_t(`difft') rho(`rho')
		
		* g_i
		tempvar g_i
		by `ivar': gen double `g_i'=1 if _n==1
		by `ivar': replace `g_i'=((1-`rho'^`difft')/(sqrt(1-`rho'^(2*`difft')))) if _n>1
		replace `g_i'=sqrt(1-`rho'^2)*`g_i'				

		*   g_i2=(g_i'g_i) 			(constant over i)
		*  g_i2s=(g_i'g_i) 			(only last obs of each panel) 
		tempvar g_i2 g_i2s 
		gen double `g_i2'=`g_i'*`g_i'
		by `ivar': gen double `g_i2s'=sum(`g_i2')
		by `ivar': replace `g_i2s'=. 			if _n<_N  
		by `ivar': replace `g_i2'=`g_i2s'[_N]

		* g_i2ss=sum_i^N(g_i'g_i) and sig_w
		tempname g_i2ss sig_w
		noi sum `g_i2s', meanonly
		scalar `g_i2ss'=r(sum)  	
		scalar `sig_w'=sqrt(`g_i2ss'*`sig_u'^2+`N_g'*`sig_e'^2)
		
		* theta_i
		if `T_con'==1 {
			tempname theta_i
			scalar `theta_i' = e(thta_50)
		}
		else {
			tempvar theta_i
			gen double `theta_i' = 1-`sig_e'/sqrt(`g_i2'*`sig_u'^2+`sig_e'^2) 
			sum `theta_i', d
			tempname theta_sum
			matrix `theta_sum'=(r(min),r(p5),r(p50),r(p95),r(max))
			matrix coln `theta_sum' = min p5 p50 p95 max
		}	
			
		* RE Baltagi-Wu transform 
		sort `ivar' `tvar'
		re_bw_transf `varlist', ivar(`ivar') thetai(`theta_i') gi(`g_i') gi2(`g_i2')
	}
}

* WALS estimates on transformed data
if "`TEMP_foc_vars'"!="" {
	#delimit;
	noi wals `depvar' (`TEMP_foc_vars') `TEMP_aux_vars'  
		`wgt_exp' if `touse'==1, 
		nocons sigma(`sigma_W') prior(`PRI') `QUADRATURE' 
		plugin(`plugin') patht(`pathtab')
		lev(`level') rep(`reps') rseed(`rseed') saving(`saving') 
		nocollcheck notable nohead `fast';
	#delimit cr
}
else {
	#delimit;
	noi wals `depvar' `TEMP_aux_vars'  
		`wgt_exp' if `touse'==1, 
		nocons sigma(`sigma_W') prior(`PRI') `QUADRATURE' 
		plugin(`plugin') patht(`pathtab')
		lev(`level') rep(`reps') rseed(`rseed') saving(`saving') 
		nocollcheck notable nohead `fast';
	#delimit cr
}	
if e(k1)!=`k1'|e(k2)!=`k2' {
	noi di as err "Transformed data may contain collinear regressors"
	exit 198
}
restore
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Parse wals results
*------------------------------------------------------------------------------------------------
tempname stddev prior_par M_W M_pm M_t_ratios
local N=e(N)
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

* Put constant term as last element  
tempname M_b 
if "`ac'"=="" {
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
		if "`ac'"=="" {
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
	if "`ac'"=="" {
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
	if "`re'"==""		{
		ereturn post `WALS_b' `wgt_exp_post', dep(`depvar') obs(`N') esample(`touse') dof(`df_r') buildfvinfo 
	}
	else {
		ereturn post `WALS_b', dep(`depvar') obs(`N') esample(`touse') buildfvinfo 
	}
	ereturn scalar n			=`N_g'
	ereturn scalar Tmin			=`g_min'
	ereturn scalar Tavg			=`g_avg'
	ereturn scalar Tmax			=`g_max'
	ereturn scalar Tcon			=`T_con'
	ereturn scalar k0 			=`LSU_k'
	ereturn scalar k1 			=`k1'
	ereturn scalar k2 			=`k2'
	if "`re'"==""	ereturn scalar df_r=`e(df_r)'			
	ereturn scalar sigma		=`stddev'
	if "`re'"!="" {
		ereturn scalar sig_nu	=`sig_u'
		ereturn scalar sig_e  	=`sig_e'
		if `e(Tcon)'==1 {
			ereturn scalar alpha=`theta_i'
		}
	}
	if "`ar1'"!=""	ereturn scalar rho=`rho'
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

	ereturn local predict 		"xtwals_p"
	ereturn local properties	"`e(properties)'"
	if "`PRI'"=="laplace" 		ereturn local quadmethod	"analytic"
	else 						ereturn local quadmethod	"`quadm'"
	ereturn local prior 		"`PRI'"
	if "`re'"=="" {
		ereturn local wexp			"`e(wexp)'"
		ereturn local wtype			"`e(wtype)'"
	}
	if "`ac'"=="" 				ereturn local constype "focus"
	else 						ereturn local constype "auxiliary"
	ereturn local auxvars 		"`M_aux_names'"
	ereturn local focvars 		"`M_foc_names'"
	ereturn local omitvars	 	"`OMITTED_vars'"
	ereturn local allvars 		"`LSU_names'"
	ereturn local depvar		"`e(depvar)'"
	if "`ar1'"!="" ereturn local tvar		"`tvar'"
	ereturn local ivar			"`ivar'"
	if "`re'"=="" & "`ar1'"=="" {
		ereturn local errors	"iid"
		ereturn local approach 	"fe"
		ereturn local model 	"fe-iid"
		ereturn local title 	"FE-WALS estimates"
	}
	else if "`re'"=="" & "`ar1'"!="" {
		ereturn local rhotype	"`rho_type'"	
		ereturn local errors	"ar1"
		ereturn local approach 	"fe"
		ereturn local model 	"fe-ar1"
		ereturn local title 	"FE-WALS estimates with AR(1) errors"
	}
	else if "`re'"!="" & "`ar1'"=="" {
		ereturn local errors	"iid"
		ereturn local approach 	"re"
		ereturn local model 	"re-iid"
		ereturn local title 	"RE-WALS estimates"
	}
	else if "`re'"!="" & "`ar1'"!="" {
		ereturn local rhotype	"`rho_type'"	
		ereturn local errors	"ar1"
		ereturn local approach 	"re"
		ereturn local model 	"re-ar1"
		ereturn local title 	"RE-WALS estimates with AR(1) errors"
	}
	ereturn local cmd 			"xtwals"

	ereturn matrix priorpar		=`prior_par'
 	ereturn matrix wals_wgt		=`M_W'
 	ereturn matrix wals_pm		=`M_pm'
 	ereturn matrix t_ratios		=`M_t_ratios'
	if "`re'"!="" & `e(Tcon)'!=1 {
		ereturn matrix alpha_stat=`theta_sum'
	}
}
else{
	if "`re'"==""		{
		ereturn post `WALS_b' `WALS_V' `wgt_exp_post', dep(`depvar') obs(`N') esample(`touse') dof(`df_r') buildfvinfo
	}
	else {
		ereturn post `WALS_b' `WALS_V', dep(`depvar') obs(`N') esample(`touse') buildfvinfo
	}
	ereturn scalar n			=`N_g'
	ereturn scalar Tmin			=`g_min'
	ereturn scalar Tavg			=`g_avg'
	ereturn scalar Tmax			=`g_max'
	ereturn scalar Tcon			=`T_con'
	ereturn scalar k0 			=`LSU_k'
	ereturn scalar k1 			=`k1'
	ereturn scalar k2 			=`k2'
	ereturn scalar rank 		=`k'
	if "`re'"==""	ereturn scalar df_r=`e(df_r)'			
	ereturn scalar sigma		=`stddev'
	if "`re'"!="" {
		ereturn scalar sig_nu	=`sig_u'
		ereturn scalar sig_e  	=`sig_e'
		if `e(Tcon)'==1 {
			ereturn scalar alpha	=`theta_i'
		}
	}
	if "`ar1'"!=""	ereturn scalar rho=`rho'
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
	ereturn scalar reps 		=`reps'
	ereturn scalar level    	=`level'

 	ereturn local marginsok 	"XB default"
	ereturn local predict 		"xtwals_p"
	ereturn local properties	"`e(properties)'"
	ereturn local bcsimdata 	`"`wals_bc_fname'"'	
	ereturn local plugin		"`BIAS_type'"
	if "`PRI'"=="laplace" 		ereturn local quadmethod	"analytic"
	else 						ereturn local quadmethod	"`quadm'"
	ereturn local prior 		"`PRI'"
	if "`re'"=="" {
		ereturn local wexp			"`e(wexp)'"
		ereturn local wtype			"`e(wtype)'"
	}
	if "`ac'"=="" 				ereturn local constype "focus"
	else 						ereturn local constype "auxiliary"
	ereturn local auxvars 		"`M_aux_names'"
	ereturn local focvars 		"`M_foc_names'"
	ereturn local omitvars	 	"`OMITTED_vars'"
	ereturn local allvars 		"`LSU_names'"
	ereturn local depvar		"`e(depvar)'"
	if "`ar1'"!="" ereturn local tvar		"`tvar'"
	ereturn local ivar			"`ivar'"
	if "`re'"=="" & "`ar1'"=="" {
		ereturn local errors	"iid"
		ereturn local approach 	"fe"
		ereturn local model 	"fe-iid"
		ereturn local title 	"FE-WALS estimates"
	}
	else if "`re'"=="" & "`ar1'"!="" {
		ereturn local rhotype	"`rho_type'"	
		ereturn local errors	"ar1"
		ereturn local approach 	"fe"
		ereturn local model 	"fe-ar1"
		ereturn local title 	"FE-WALS estimates with AR(1) errors"
	}
	else if "`re'"!="" & "`ar1'"=="" {
		ereturn local errors	"iid"
		ereturn local approach 	"re"
		ereturn local model 	"re-iid"
		ereturn local title 	"RE-WALS estimates"
	}
	else if "`re'"!="" & "`ar1'"!="" {
		ereturn local rhotype	"`rho_type'"	
		ereturn local errors	"ar1"
		ereturn local approach 	"re"
		ereturn local model 	"re-ar1"
		ereturn local title 	"RE-WALS estimates with AR(1) errors"
	}
	ereturn local cmd 			"xtwals"
	
	ereturn matrix priorpar		=`prior_par'
 	ereturn matrix wals_wgt		=`M_W'
 	ereturn matrix wals_pm		=`M_pm'
 	ereturn matrix t_ratios		=`M_t_ratios'
	if "`re'"!="" & `e(Tcon)'!=1 {
		ereturn matrix alpha_stat=`theta_sum'
	}
	ereturn matrix ci			=`WALS_CI'
	ereturn matrix MSE			=`WALS_MSE'
	ereturn matrix rmse			=`WALS_rmse'
	ereturn matrix bias			=`WALS_bias'
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Display estimation results
*------------------------------------------------------------------------------------------------
if "`fast'"=="" 	noi Di_xtwals, level(`level') `display_options'   
*------------------------------------------------------------------------------------------------

}
end
*-----------------------------------------------------------------------------------------------------------------------	



*-----------------------------------------------------------------------------------------------------------------------	
* FE transform
*-----------------------------------------------------------------------------------------------------------------------	
program define fe_transf
	syntax varlist [aw fw], ivar(varname)
	tokenize `varlist'
	while ("`1'"!="") {
		tempvar tv
		by `ivar': gen double `tv' = sum(`1')/_n
		summ `1' [`weight'`exp']
		by `ivar': replace `tv' = `1' - `tv'[_N] + r(mean)
		drop `1'
		rename `tv' `1'
		mac shift
	}
end	
*-----------------------------------------------------------------------------------------------------------------------	



*-----------------------------------------------------------------------------------------------------------------------	
* RE transform
*-----------------------------------------------------------------------------------------------------------------------	
program define re_transf
	syntax varlist, ivar(varname) thetai(string)
	tempvar tv
	tokenize `varlist'
	while ("`1'"!="") {
		by `ivar': gen double `tv' = sum(`1')/_n
		by `ivar': replace `tv' = `1' - `thetai'*`tv'[_N] 
		drop `1'
		rename `tv' `1'
		mac shift
	}
end	
*-----------------------------------------------------------------------------------------------------------------------	



*-----------------------------------------------------------------------------------------------------------------------	
* RE Baltagi-Wu transform
*-----------------------------------------------------------------------------------------------------------------------	
program define re_bw_transf
	syntax varlist, ivar(varname) thetai(string) gi(varname) gi2(varname)
	tempvar tv 
	tokenize `varlist'
	while ("`1'"!="") {
		gen double `tv'=`1'*`gi'
		by `ivar': replace `tv'=sum(`tv')
		by `ivar': replace `tv'=`tv'[_N]
		replace `1'=`1'-`thetai'*`gi'*(`tv'/`gi2')
		capture drop `tv'
		mac shift
	}	
end	
*-----------------------------------------------------------------------------------------------------------------------	



*-----------------------------------------------------------------------------------------------------------------------	
* C_i(rho) transform 
*-----------------------------------------------------------------------------------------------------------------------	
program define Ci_rho_transf
	syntax varlist, ivar(varname) dif_t(varname) rho(string)
	tempvar tv
	tokenize `varlist'
	while ("`1'"!="") {
		by `ivar': gen double `tv'=(sqrt(1-`rho'^2))*`1' 	if _n==1
		by `ivar': replace `tv'=(sqrt(1-`rho'^2))*			///
			(`1'[_n]*(1/sqrt((1-`rho'^(2*`dif_t')))) - 		///
			`1'[_n-1]*(`rho'^(`dif_t')/ 					///
				sqrt(1-`rho'^(2*`dif_t')))) 				if _n>1
		drop `1'
		rename `tv' `1'
		mac shift
	}
end	
*-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// Di_xtwals: display of xtwals output 
//-----------------------------------------------------------------------------------------------------------------------	
program Di_xtwals
	syntax [, LEVel(cilevel) noTABle noFOCus noAUXiliary noHEADer cformat(string)]

	local y=abbrev(`"`e(depvar)'"',10)
	local focvars 	"`e(focvars)'"
	local auxvars	"`e(auxvars)'"
	local k1=e(k1)
	local k2=e(k2)
	local regressors: coln 	 e(b)
	local LSU_k		: colsof e(b)
	local ERR "`e(errors)'"
	local cons "`e(constype)'"
	if "`cons'"=="focus" 	{
		local k1_nc=`k1'-1
		local k2_nc=`k2'
		local focvars_nc: subinstr local focvars "_cons" ""
		local M_names "`focvars_nc' `auxvars' _cons"
	}
	else {
		local k1_nc=`k1'
		local k2_nc=`k2'-1
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
		if "`ERR'"=="iid" {
			#delimit;
			di as text _n `"`e(title)'"'					
				   as text _col(`p1') "Number of obs" 	_col(`p2') "="
				   as res %8.0f e(N)
				
				_n as text "Prior            : " 
				   as res "`e(prior)'" 
				   as text _col(`p1') "k1 " 			_col(`p2') "=" 	
				   as res %8.0f e(k1) 							

				_n as text "Panel variable   : " 	 
				   as res e(ivar) 							
				   as text _col(`p1') "k2 " 			_col(`p2') "=" 	
				   as res %8.0f e(k2) 									

				_n as text "Number of panels : " 	 
				   as res e(n) 							
				   as text _col(`p1') "sigma " 			_col(`p2') "=" 
				   as res %8.0g e(sigma) 
				   
				_n as text "Obs per panel    : min=" 
				   as res `e(Tmin)'	 
				   as text ", max=" 
				   as res `e(Tmax)'	
				   as text _col(`p1') "MC reps " 		_col(`p2') "=" 
				as res %8.0f e(reps); 								
			#delimit cr	
		}
		else {
			#delimit;
			di as text _n `"`e(title)'"'					
				   as text _col(`p1') "Number of obs" 	_col(`p2') "="
				   as res %8.0f e(N)
				
				_n as text "Prior            : " 
				   as res "`e(prior)'" 
				   as text _col(`p1') "k1 " 			_col(`p2') "=" 	
				   as res %8.0f e(k1) 				

				_n as text "Panel variable   : " 	 
				   as res e(ivar) 							
				   as text _col(`p1') "k2 " 			_col(`p2') "=" 	
				   as res %8.0f e(k2) 									
				   
				_n as text "Time variable    : " 	 
				   as res e(tvar) 				
				   as text _col(`p1') "rho " 			_col(`p2') "=" 
				   as res %8.0g e(rho) 

				_n as text "Number of panels : "  
				   as res e(n) 				
				   as text _col(`p1') "sigma " 			_col(`p2') "=" 
				   as res %8.0g e(sigma) 
				   
				_n as text "Obs per panel    : min=" 
				   as res `e(Tmin)'	 
				   as text ", max=" 
				   as res `e(Tmax)'	
				   as text _col(`p1') "MC reps " 		_col(`p2') "=" 
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
			if `r(omit)'==0 & _b[`xx']!=0 {
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

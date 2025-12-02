* Weighted Average Least Squares (WALS) for GLM
*! Version 1.0
* Date			: 2024/10/05
* Authors		: De Luca Giuseppe, Jan R. Magnus
*-------------------------------------------------------------------------------------------------
* glmwals
*-------------------------------------------------------------------------------------------------
program define glmwals, eclass sort
	version 14.0, missing
	if replay() {
		if "`e(cmd)'" != "glmwals" {
			error 301
		}
		syntax [, SAVing(string asis) noHEADer noTABle noFOCus noAUXiliary LEVel(cilevel) cformat(string) eform]

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
		Di_glmwals, `header' `table' `focus' `auxiliary' lev(`level') `cformat' `eform'
		exit
	}
	if !replay() {
		syntax [anything] [fw aw iw /] [if] [in], [*]
		local cmdline : copy local 0
	}

	Estimate `0'
	ereturn local cmdline `"glmwals `0'"'
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
cap macro drop GLMWALS*
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
	sigma(real 0)  											

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
	
	Family(string)
	Link(string)
	OFFset(varname numeric)			
	LNOFFset(varname numeric)		
	EXPosure(varname numeric)		
	INIt(string)
	LTOLerance(real 1e-6)
	iter0(integer 300) 
	onestep	
	TOLerance(real 1e-6)
	ITERate(integer 300) 
	noLOG
	eform 					
	SHOWIrls
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
* Mark estimation sample 
*------------------------------------------------------------------------------------------------
marksample touse
markout `touse' `depvar' `foclist0' `auxlist0' `exp' 
count if `touse' 
if r(N)==0 	error 2000 
if r(N)==1 	error 2001
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Preliminary check of options
*------------------------------------------------------------------------------------------------
global GLMWALS_y	`"`depvar'"'		

* init
GLMWALS_Map_mu0 `"`init'"'
local start `r(start)'
if "`start'"=="3" 	{
	local mu0_var `r(mu0_var)'
	local init "init(`mu0_var')"
}	
else {				
	local mu0_var
	local init
}

* ltolerance
if `ltolerance'>=1 | `ltolerance'<=0 {
	di as err "{bf:ltolerance} must be in (0,1)"
	exit 198
}

* iter0
if `iter0' < 1 {
	di as err "{bf:iter0} must be a positive integer"
	exit 198
}

* iterate
if `iterate' < 1 {
	di as err "{bf:iterate} must be a positive integer"
	exit 198
}

* tolerance
if `tolerance'>=1 | `tolerance'<=0 {
	di as err "{bf:tolerance} must be in (0,1)"
	exit 198
}

* nolog
if "`log'"=="nolog" local log_yn "`log'"
else			 	local log_yn "log"

* offset
if "`offset'" != "" & "`lnoffset'" != "" {
	di as err "cannot specify both {bf:offset} and {bf:lnoffset}"
	exit 198
}
if "`offset'" != "" & "`exposure'" != "" {
	di as err "cannot specify both {bf:offset} and {bf:exposure}"
	exit 198
}
if "`lnoffset'" != "" & "`exposure'" != "" {
	di as err "cannot specify both {bf:lnoffset} and {bf:exposure}"
	exit 198
}
if "`exposure'" != "" {
	local lnoffset `exposure'
	local exposure
}
if "`lnoffset'" != "" {
	capture assert `lnoffset' > 0 if `touse'
	if _rc {
		di as err "{bf:exposure} must be greater than zero"
		exit 459
	}
	tempvar offset
	qui gen double `offset' = ln(`lnoffset')
	local offvar "ln(`lnoffset')"
}
if "`offset'" != "" {
	markout `touse' `offset'
	local offopt "offset(`offset')"
	if "`offvar'" == "" {
		local offvar "`offset'"
	}
	local moffset = "-`offset'"
}

* glm options 
#delimit;
	local glm_opts_rm "irls vce(eim) `offopt' 
		`init' iterate(`iter0') ltolerance(`ltolerance')  
		nolog noheader notable";
#delimit cr
if "`link'"  != ""  		local glm_opts_rm `"`glm_opts_rm' link(`link')"'
if "`family'"!= "" 			local glm_opts_rm `"`glm_opts_rm' family(`family')"'
if "`showirls'"=="" {
	local glm_opts `glm_opts_rm'
}
else {
	#delimit;
		local glm_opts "irls vce(eim) `offopt' 
			`init' iterate(`iter0') ltolerance(`ltolerance')
			`log_yn' `eform'";
	#delimit cr
	if "`link'"  != ""  		local glm_opts `"`glm_opts' link(`link')"'
	if "`family'"!= "" 			local glm_opts `"`glm_opts' family(`family')"'
}

* onestep  
if "`onestep'"!="" & `tolerance'!=1e-6 {
	noi di as err "cannot specify both {bf:onestep} and {bf:tolerance}"
	exit 198
}
if "`onestep'"!="" & `iterate'!=300 {
	noi di as err "cannot specify both {bf:onestep} and {bf:iterate}"
	exit 198
}
if "`onestep'"!="" & "`log_yn'"=="nolog" {
	noi di as err "cannot specify both {bf:onestep} and {bf:nolog}"
	exit 198
}

* Display options
if "`fast'"=="" {
	if "`cformat'"!="" local cformat "cformat(`cformat')"
	local display_options "`table' `focus' `auxiliary' `header' `cformat' `eform'"
}

* saving
cap erase `"`c(tmpdir)'/STWALS_000001.tmp"'
noi _savingopt_parse CIMC_save CIMC_savereplace : saving ".dta" `"`saving'"'
if `"`CIMC_save'"'!="" & `"`CIMC_savereplace'"'=="" confirm new file `"`CIMC_save'"'
if `"`CIMC_save'"'=="" {
	local CIMC_save `"`c(tmpdir)'/STWALS_000001.tmp"'
	local CIMC_savereplace "replace"
}
*------------------------------------------------------------------------------------------------


	
*------------------------------------------------------------------------------------------------
* Weights 
*------------------------------------------------------------------------------------------------
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
	gen byte `wgt_var' = `touse' 			if `touse'
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
* note: if quadext was not specified, this is created here to speed up the iterative algorithm 
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
	else {
		tempname wals_GL_PW
		if "`PRI'"=="horseshoe" 						mata: `wals_GL_PW'=gausslegendre(`quadnpts')
		if "`PRI'"!="laplace" & "`PRI'"!="horseshoe" 	mata: `wals_GL_PW'=gausslaguerre(`quadnpts')
		local quadext "`wals_GL_PW'"
	}
	if `quadatol'!=1e-9   	noi di as text "note: {bf:quadatol} is ineffective with {bf:quadmethod(gauss)}"
	if `quadrtol'!=1e-7		noi di as text "note: {bf:quadrtol} is ineffective with {bf:quadmethod(gauss)}"
	if  "`PRI'"!="laplace" 	local QUADRATURE "quadm(`quadmethod') quadext(`quadext')"
	else 					local QUADRATURE "quadm(`quadmethod')"
}
if "`quadm'"=="adaptive" {		
	if "`quadext'"!=""		noi di as text "note: {bf:quadext} is ineffective with {bf:quadmethod(adaptive)}"
	if "`quadnpts'"!="500"	noi di as text "note: {bf:quadnpts} is ineffective with {bf:quadmethod(adaptive)}"
	if  "`PRI'"!="laplace" 	local QUADRATURE "quadm(`quadmethod') quadat(`quadatol') quadrt(`quadrtol')"
	else 					local QUADRATURE "quadm(`quadmethod')"
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
* Family and link																			TBA
*------------------------------------------------------------------------------------------------
GLMWALS_MapFL 	`"`family'"' `"`link'"'
local family    `"`r(family)'"'
local link      `"`r(link)'"'
local pow       `r(power)'
local m         `r(m)'
local mfixed    `r(mfixed)'
if "`family'" == ""|"`link'" == "" {
	di as err "incomplete specification of family() and link()"
	exit 198
}
if "`family'"=="glmwals_v1" & "`link'"=="glmwals_l01" local onestep "onestep"

if "`m'" != "" {
	capture confirm number `m'
	if _rc!=0 	markout `touse' `m' 
}
global GLMWALS_V	`"`family'"'		/* V(mu) program 		*/
global GLMWALS_L	`"`link'"'			/* g(mu) program 		*/
global GLMWALS_m	`"`m'"'				/* Binomial denominator */
global GLMWALS_p	`"`pow'"'			/* Power   			 	*/
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Check for perfect collinearity in the unrestricted model
*------------------------------------------------------------------------------------------------
* IRLS estimates of unrestricted model	
if "`showirls'"!="" noi di as text "IRLS estimates of unrestricted model"
cap noi glm `depvar' `foclist0' `auxlist0' `wgt_exp' if `touse', `nc' `glm_opts' 
if _rc!=0 {
	noi di as err "IRLS estimates of unrestricted model do not converge"
	exit `e(rc)'
}
local MLU_k: colsof e(b)
local MLU_names: coln e(b)

* Estimation sample from unrestricted model
replace `touse'=e(sample)
sum `touse' `wgt_exp' if `touse', meanonly
local N=int(r(sum))

* List of noncollinear regressors, excluding the constant term
_ms_omit_info e(b)
tempname x0_omit
matrix `x0_omit'=r(omit)
local x_names ""
local jj=1
foreach xx of local MLU_names {
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
if "`nc'"~="noconstant" & "`ac'"==""	local TEMP_foc_vars `const' `TEMP_foc_vars' 
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
	noi di as err "Model must contain at least one auxiliary regressor (including, if any, the constant term)" 
	error 102
}
local k = `k1'+`k2'
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Scale parameter 
*------------------------------------------------------------------------------------------------
if `sigma'<0 {
	di as err "{bf:sigma} must be nonnegative"
	error 411
}
if `sigma'==0 {
	if `"`family'"'=="glmwals_v2"|`"`family'"'=="glmwals_v3" {
		local sigma_MLU=1
	}
	else {
		local sigma_MLU=e(dispers_p)^.5
	}
	cap assert `sigma_MLU'!=.
	if _rc {
		noi di as err "Pearson dispersion parameter is missing"
		error 9
	}
}
else {
	local sigma_MLU=`sigma'
	cap assert `sigma_MLU'!=.
	if _rc {
		noi di as err "{bf:sigma} cannot be missing"
		error 9
	}
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Starting values 
*------------------------------------------------------------------------------------------------
* Set titles and check range of depvar (and m) 
$GLMWALS_V -1 `touse'
$GLMWALS_L -1 `touse'

* Initial values for bar_mu and bar_eta
tempvar bar_mu bar_eta
if "`start'"=="1" {
	predict double `bar_mu'														if `touse', mu
	$GLMWALS_mu `bar_mu'	
	cap assert `bar_mu'!=. 														if `touse'
	if _rc!=0 {
		noi di as err "Missing values in the initial values for mu"
		error 9
	}
	predict double `bar_eta'													if `touse', xb		
	cap assert `bar_eta'!=. 													if `touse'
	if _rc!=0 {
		noi di as err "Missing values in the initial values for eta"
		error 9
	}
}
if "`start'"=="2" {
	if "`nc'"~="noconstant" & "`ac'"=="" {
		cap glm `depvar' `foclist0' `wgt_exp' if `touse', `glm_opts_rm' 
	}
	else {	
		cap glm `depvar' `foclist0' `wgt_exp' if `touse', `glm_opts_rm' noconstant   
	}
	if _rc!=0 {
		noi di as err "IRLS estimates of restricted model do not converge"
		exit `e(rc)'
	}
	predict double `bar_mu'														if `touse', mu
	$GLMWALS_mu `bar_mu'	
	cap assert `bar_mu'!=. 														if `touse'
	if _rc!=0 {
		noi di as err "Missing values in the initial values for mu"
		error 9
	}
	predict double `bar_eta'													if `touse', xb		
	cap assert `bar_eta'!=. 													if `touse'
	if _rc!=0 {
		noi di as err "Missing values in the initial values for eta"
		error 9
	}
}
if "`start'"=="3" {
	if `"`mu0_var'"'!="" {
		gen double `bar_mu' = `mu0_var' 										if `touse'
	}
	else {
		sum `depvar' `wgt_exp' if `touse', mean
		if "${GLMWALS_V}"!="glmwals_v2" {
			gen double `bar_mu' = (`depvar'+r(mean))/(`m'+1)
		}
		else {
			gen double `bar_mu' = `m'*(`depvar'+.5)/(`m'+1)
		}
	}
	$GLMWALS_mu `bar_mu'	
	cap assert `bar_mu'!=. 														if `touse'
	if _rc!=0 {
		noi di as err "Missing values in the initial values for mu"
		error 9
	}
	$GLMWALS_L 0 `bar_eta' `bar_mu'												// it does not include offset
	cap assert `bar_eta'!=. 													if `touse'
	if _rc!=0 {
		noi di as err "Missing values in the initial values for eta"
		error 9
	}
}

* Initial values for the derivative of the inverse link function and the variance function  
tempvar bar_dmu bar_v 
$GLMWALS_L 2 `bar_eta' `bar_mu' `bar_dmu'										/* (d mu)/(d eta) 	*/
cap assert `bar_dmu'!=. 														if `touse'
if _rc!=0 {
	noi di as err "Missing values in the initial values for (d mu)/(d eta)"
	error 9
}
$GLMWALS_V 1 `bar_eta' `bar_mu' `bar_v'											/* v(mu) 			*/
cap assert `bar_v'!=. 															if `touse'
if _rc!=0 {
	noi di as err "Missing values in the initial values for v(mu)"
	error 9
}

* Initial values for the weights
tempvar bar_w bar_w_n 
gen double `bar_w' = `bar_dmu'*`bar_dmu'/`bar_v' 								if `touse'
cap assert `bar_w'!=. 															if `touse'
if _rc!=0 {
	noi di as err "Missing values in the initial values for the weights"
	error 9
}
cap assert `bar_w'>=0 															if `touse'
if _rc!=0 {
	noi di as err "Negative values in the initial values for the weights"
	error 411
}
sum `bar_w' `wgt_exp', meanonly
local bar_w_m=r(mean)
gen double `bar_w_n' = `wgt_var' *`bar_w'/`bar_w_m' 							if `touse'

* Initial values for sigma in wals  
local sigma_W=`sigma_MLU'/`bar_w_m'^.5

* Initial values for adjusted depvar
tempvar bar_depvar
if "`start'"=="1"|"`start'"=="2" {
	gen double `bar_depvar'=`bar_eta' `moffset' + (`depvar'-`bar_mu')/`bar_dmu' if `touse'
}
else{
	gen double `bar_depvar'=`bar_eta' + (`depvar'-`bar_mu')/`bar_dmu' 			if `touse'
}
cap assert `bar_depvar'!=. 														if `touse'
if _rc!=0 {
	noi di as err "Missing values in the initial values for the adjusted depvar"
	error 9
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* WALS estimates 
*------------------------------------------------------------------------------------------------
if "`onestep'"!="" {
	* Title
	local glmwals_tit "One-step WALS estimates of GLM"

	* One-step WALS estimates (all calculations)
	if "`TEMP_foc_vars'"!="" {
		#delimit;
		wals `bar_depvar' (`TEMP_foc_vars') `TEMP_aux_vars' if `touse'==1 
			[iw=`bar_w_n'], 
			nocons sigma(`sigma_W') prior(`PRI') `QUADRATURE' 
			plugin(`plugin') patht(`pathtab')
			lev(`level') rep(`reps') rseed(`rseed') saving(`saving') 
			nocollcheck notable nohead `fast';
		#delimit cr
	}
	else {
		#delimit;
		wals `bar_depvar' `TEMP_aux_vars' if `touse'==1 
			[iw=`bar_w_n'], 
			nocons sigma(`sigma_W') prior(`PRI') `QUADRATURE' 
			plugin(`plugin') patht(`pathtab')
			lev(`level') rep(`reps') rseed(`rseed') saving(`saving') 
			nocollcheck notable nohead `fast';
		#delimit cr
	}
	if e(k1)!=`k1'|e(k2)!=`k2' {
		noi di as err "Collinear regressors in the transformed data"
		exit 9
	}
	local converged =.
	local glmwals_iter=.   
}
else {
	* Title
	local glmwals_tit "Iterative WALS estimates of GLM"
	if "`showirls'"!="" noi di _n as text "`glmwals_tit'"
	
	* One-step WALS estimates (point estimates only)
	if "`TEMP_foc_vars'"!="" {
		#delimit;
		wals `bar_depvar' (`TEMP_foc_vars') `TEMP_aux_vars' if `touse'==1 
			[iw=`bar_w_n'], 
			nocons sigma(`sigma_W') prior(`PRI') `QUADRATURE' 
			nocollcheck fast;
		#delimit cr
	}
	else {
		#delimit;
		wals `bar_depvar' `TEMP_aux_vars' if `touse'==1 
			[iw=`bar_w_n'], 
			nocons sigma(`sigma_W') prior(`PRI') `QUADRATURE' 
			nocollcheck fast;
		#delimit cr
	}
	if e(k1)!=`k1'|e(k2)!=`k2' {
		noi di as err "Collinear regressors in the transformed data"
		exit 9
	}
	if "`log_yn'"!="nolog" 	noi di _n as txt "Iteration 1:  Rel. Diff. e(b) initialized" 

	* Iterative WALS estimates (point estimates only)
	tempname bar_b new_b mrdif_b
	scalar `mrdif_b'=99
	local space " "
	local glmwals_iter=1   
	while (`mrdif_b'>`tolerance' & `glmwals_iter'<`iterate') {
		
		* Drop data transformations from previous step
		cap drop `bar_eta' `bar_mu' `bar_dmu' `bar_v' `bar_w' `bar_w_n' `bar_depvar'  

		* Update estimated coeffcients
		matrix `bar_b'=e(b)

		* Update bar_eta 
		gen double `bar_eta'=0													if `touse'
		local ss=1
		foreach xx of local TEMP_foc_vars {
			replace `bar_eta'=`bar_eta' + `bar_b'[1,`ss'] * `xx'				if `touse'
			local ss=`ss'+1
		}
		foreach xx of local TEMP_aux_vars {
			replace `bar_eta'=`bar_eta' + `bar_b'[1,`ss'] * `xx'				if `touse'
			local ss=`ss'+1
		}
		cap assert `bar_eta'!=. 												if `touse'
		if _rc!=0 {
			noi di as err "Missing values in eta (itearation `glmwals_iter')"
			error 9
		}
		if "`offset'" != "" {
			replace `bar_eta' = `bar_eta'+`offset'
		}

		* Update bar_mu 
		$GLMWALS_L 1 `bar_eta' `bar_mu' 							/* mu = g^-1(eta) 	*/
		$GLMWALS_mu `bar_mu'	
		cap assert `bar_mu'!=. 													if `touse'
		if _rc!=0 {
			noi di as err "Missing values in mu (itearation `glmwals_iter')"
			error 9
		}

		* Update bar_dmu and bar_v   
		$GLMWALS_L 2 `bar_eta' `bar_mu' `bar_dmu'					/* (d mu)/(d eta) 	*/
		$GLMWALS_V 1 `bar_eta' `bar_mu' `bar_v'						/* v(mu) 			*/
		
		* Update bar_w
		gen double `bar_w' = `bar_dmu'*`bar_dmu'/`bar_v' 						if `touse'
		cap assert `bar_w'!=. 													if `touse'
		if _rc!=0 {
			noi di as err "Missing weights (iteration `glmwals_iter')"
			error 9
		}
		cap assert `bar_w'>=0 													if `touse'
		if _rc!=0 {
			noi di as err "Negative weights (iteration `glmwals_iter')"
			error 411
		}
		sum `bar_w' `wgt_exp', meanonly
		local bar_w_m=r(mean)
		gen double `bar_w_n' = `wgt_var'*`bar_w'/`bar_w_m' 						if `touse'

		* Unpdate sigma in wals  
		local sigma_W=`sigma_MLU'/`bar_w_m'^.5
		
		* Update bar_depvar
		gen double `bar_depvar'=`bar_eta' `moffset' + (`depvar'-`bar_mu')/`bar_dmu' if `touse'
		cap assert `bar_depvar'!=. 													if `touse'
		if _rc!=0 {
			noi di as err "Missing values in adjusted depvar (itearation `glmwals_iter')"
			error 9
		}
		
		* New WALS estimates (point estimates only)  
		if "`TEMP_foc_vars'"!="" {
			#delimit;
			noi wals `bar_depvar' (`TEMP_foc_vars') `TEMP_aux_vars' if `touse'==1 
				[iw=`bar_w_n'], 
				nocons sigma(`sigma_W') prior(`PRI') `QUADRATURE' nocollcheck fast;
			#delimit cr
		}
		else {
			#delimit;
			noi wals `bar_depvar' `TEMP_aux_vars' if `touse'==1 
				[iw=`bar_w_n'], 
				nocons sigma(`sigma_W') prior(`PRI') `QUADRATURE' nocollcheck fast;
			#delimit cr
		}
		if e(k1)!=`k1'|e(k2)!=`k2' {
		noi di as err "Collinear regressors in the transformed data"
			exit 9
		}
		matrix `new_b'=e(b)

		* Update iteration
		local glmwals_iter=`glmwals_iter'+1
		if `glmwals_iter'>=10 local space ""
		
		* mreldif e(b)
		scalar `mrdif_b'=mreldif(`bar_b',`new_b')
		if "`log_yn'"!="nolog" noi di as txt "Iteration `glmwals_iter':`space' Rel. Diff. e(b) = " as res %9.0g `mrdif_b'
	}
	if `glmwals_iter'>=`iterate'|`mrdif_b'>`tolerance' 	local converged=0
	else												local converged=1

	* Iterative WALS estimates (all calculations)
	if "`TEMP_foc_vars'"!="" {
		#delimit;
		noi wals `bar_depvar' (`TEMP_foc_vars') `TEMP_aux_vars' if `touse'==1 
			[iw=`bar_w_n'],  
			nocons sigma(`sigma_W') prior(`PRI') `QUADRATURE' 
			plugin(`plugin') patht(`pathtab')
			lev(`level') rep(`reps') rseed(`rseed') saving(`saving')  
			nocollcheck notable nohead `fast';
		#delimit cr
	}
	else {
		#delimit;
		noi wals `bar_depvar' `TEMP_aux_vars' if `touse'==1 
			[iw=`bar_w_n'],  
			nocons sigma(`sigma_W') prior(`PRI') `QUADRATURE' 
			plugin(`plugin') patht(`pathtab')
			lev(`level') rep(`reps') rseed(`rseed') saving(`saving')  
			nocollcheck notable nohead `fast';
		#delimit cr
	}
	if e(k1)!=`k1'|e(k2)!=`k2' {
		noi di as err "Collinear regressors in the transformed data"
		exit 9
	}
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Parse wals results
*------------------------------------------------------------------------------------------------
tempname prior_par M_W M_pm M_t_ratios
matrix `prior_par'=e(priorpar)
matrix `M_W'=e(wals_wgt)
matrix `M_pm'=e(wals_pm)
matrix `M_t_ratios'=e(t_ratios)
if "`fast'"=="" {
	local level=e(level)
	local reps=e(reps)
	local BIAS_type	"`e(plugin)'"
}
local df_r=`N'-`k'
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
matrix `WALS_b'=J(1,`MLU_k',0)
local ss=1
forvalues jj=1(1)`MLU_k'{
	local xx: word `jj' of `MLU_names'
	local yy: word `ss' of `M_names'
	if "`xx'"=="`yy'" {
		matrix `WALS_b'[1,`jj']=`M_b'[1,`ss']
		local ++ss
	}
}
matrix coln `WALS_b'=`MLU_names'
if "`fast'"=="" {
	tempname WALS_bias WALS_V WALS_MSE WALS_rmse

	foreach mm in bias rmse {
		matrix `WALS_`mm''=J(1,`MLU_k',0)
		local ss=1
		forvalues jj=1(1)`MLU_k'{
			local xx: word `jj' of `MLU_names'
			local yy: word `ss' of `M_names'
			if "`xx'"=="`yy'" {
				matrix `WALS_`mm''[1,`jj']=`M_`mm''[1,`ss']
				local ++ss
			}
		}
		matrix coln `WALS_`mm''=`MLU_names'
	}
	foreach mm in V MSE {
		matrix `WALS_`mm''=J(`MLU_k',`MLU_k',0)
		local s_rr=1
		forvalues rr=1(1)`MLU_k'{
			local xx_rr: word `rr'   of `MLU_names'
			local yy_rr: word `s_rr' of `M_names'
			local s_jj=1
			forvalues jj=1(1)`MLU_k'{
				local xx_jj: word `jj'   of `MLU_names'
				local yy_jj: word `s_jj' of `M_names'
				if "`xx_rr'"=="`yy_rr'" & "`xx_jj'"=="`yy_jj'" {
					matrix `WALS_`mm''[`rr',`jj']=`M_`mm''[`s_rr',`s_jj']
					local ++s_jj
				}
			}
			if "`xx_rr'"=="`yy_rr'" local ++s_rr
		}
		matrix coln `WALS_`mm''=`MLU_names'
		matrix rown `WALS_`mm''=`MLU_names'
	}
}

* Initialize confidence intervals
if "`fast'"=="" {
	tempname WALS_CI	
	matrix `WALS_CI'=J(2,`MLU_k',.)
	matrix rown `WALS_CI'= ll ul
	matrix coln `WALS_CI'=`MLU_names'
}

* List of omitted regressors
local M_focaux_names `M_foc_names' `M_aux_names' 
local OMITTED_vars: list MLU_names - M_focaux_names
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
	restore
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Return estimation results
*------------------------------------------------------------------------------------------------
if "`fast'"!="" {
	ereturn post `WALS_b' `wgt_exp_post', dep(`depvar') obs(`N') esample(`touse') dof(`df_r') buildfvinfo
	ereturn scalar k0 			=`MLU_k'
	ereturn scalar k1 			=`k1'
	ereturn scalar k2 			=`k2'
	ereturn scalar sigma		=`sigma_MLU'
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
	ereturn scalar ic 			=`glmwals_iter'
	ereturn scalar converged 	=`converged'
	ereturn scalar power 		=$GLMWALS_p	
	
	ereturn local predict 		"glmwals_p"
	ereturn local properties	"`e(properties)'"
	if "`PRI'"=="laplace" 		ereturn local quadmethod	"analytic"
	else 						ereturn local quadmethod	"`quadm'"
	ereturn local prior 		"`PRI'"
	ereturn local wexp			"`e(wexp)'"
	ereturn local wtype			"`e(wtype)'"
	ereturn local offset  		"`offvar'"
	if "`nc'"=="noconstant"						ereturn local constype "`nc'"
	else if "`nc'"~="noconstant" & "`ac'"=="" 	ereturn local constype "focus"
	else 										ereturn local constype "auxiliary"
	ereturn local auxvars 		"`M_aux_names'"
	ereturn local focvars 		"`M_foc_names'"
	ereturn local omitvars	 	"`OMITTED_vars'"
	ereturn local allvars 		"`MLU_names'"
	if "`onestep'" =="" 		ereturn local esttype "iterative"
	else 						ereturn local esttype "one-step"
	ereturn local btrials 		"$GLMWALS_m"
	ereturn local linkf 		"$GLMWALS_lf"
	ereturn local linkt 		"$GLMWALS_lt"
	ereturn local link 			"$GLMWALS_L"
	ereturn local varfuncf 		"$GLMWALS_vf"
	ereturn local varfunct 		"$GLMWALS_vt"
	ereturn local varfunc 		"$GLMWALS_V"
	ereturn local title 		"`glmwals_tit'"
	ereturn local cmd 			"glmwals"	
	
	ereturn matrix priorpar		=`prior_par'
	ereturn matrix wals_wgt		=`M_W'
	ereturn matrix wals_pm		=`M_pm'
	ereturn matrix t_ratios		=`M_t_ratios'
}
else{
	ereturn post `WALS_b' `WALS_V' `wgt_exp_post', dep(`depvar') obs(`N') esample(`touse') dof(`df_r') buildfvinfo
	ereturn scalar k0 			=`MLU_k'
	ereturn scalar k1 			=`k1'
	ereturn scalar k2 			=`k2'
	ereturn scalar rank			=`k'
	ereturn scalar sigma		=`sigma_MLU'
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
	ereturn scalar ic 			=`glmwals_iter'
	ereturn scalar converged 	=`converged'
	ereturn scalar power 		=$GLMWALS_p

	ereturn local marginsok 	"mu default"
	ereturn local predict 		"glmwals_p"
	ereturn local properties	"`e(properties)'"
	ereturn local bcsimdata 	`"`CIMC_save'"'	
	ereturn local plugin		"`BIAS_type'"
	if "`PRI'"=="laplace" 		ereturn local quadmethod	"analytic"
	else 						ereturn local quadmethod	"`quadm'"
	ereturn local prior 		"`PRI'"
	ereturn local wexp			"`e(wexp)'"
	ereturn local wtype			"`e(wtype)'"
	ereturn local offset  		"`offvar'"
	if "`nc'"=="noconstant"						ereturn local constype "`nc'"
	else if "`nc'"~="noconstant" & "`ac'"==""	ereturn local constype "focus"
	else 										ereturn local constype "auxiliary"
	ereturn local auxvars 		"`M_aux_names'"
	ereturn local focvars 		"`M_foc_names'"
	ereturn local omitvars	 	"`OMITTED_vars'"
	ereturn local allvars 		"`MLU_names'"
	if "`onestep'" =="" 		ereturn local esttype "iterative"
	else 						ereturn local esttype "one-step"
	ereturn local btrials 		"$GLMWALS_m"
	ereturn local linkf 		"$GLMWALS_lf"
	ereturn local linkt 		"$GLMWALS_lt"
	ereturn local link 			"$GLMWALS_L"
	ereturn local varfuncf 		"$GLMWALS_vf"
	ereturn local varfunct 		"$GLMWALS_vt"
	ereturn local varfunc 		"$GLMWALS_V"
	ereturn local depvar		"`e(depvar)'"
	ereturn local title 		"`glmwals_tit'"
	ereturn local cmd 			"glmwals"	
	
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
cap macro drop GLMWALS*
if "`fast'"=="" 	noi Di_glmwals, level(`level') `display_options'   
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Check convergence of iterative WALS estimates
*------------------------------------------------------------------------------------------------
if "`onestep'"=="" & "`converged'"=="0" {
	noi di in red "Convergence not achieved"
	exit 430
}
*------------------------------------------------------------------------------------------------

}
end
*-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// GLMWALS_Map_mu0: map the method for initial values
//-----------------------------------------------------------------------------------------------------------------------	
program GLMWALS_Map_mu0, rclass 
	args init0 
	tokenize `"`init0'"'
	if `"`3'"'!="" {
		di as err "init(`init0') invalid"
		exit 198
	}
	local meth = lower(trim(`"`1'"'))
	local s1
	local mu0_var
	if `"`meth'"'=="" { 
		local s1 "1" 
	}
	else if `"`meth'"'=="irls_u" {
		local s1 "1" 
		if `"`2'"'!="" {
			di as err "init(`init0') invalid. `meth' does not require additional arguments"
			exit 198
		}
	}
	else if `"`meth'"'=="irls_r" {
		local s1 "2" 
		if `"`2'"'!="" {
			di as err "init(`init0') invalid. `meth' does not require additional arguments"
			exit 198
		}
	}
	else if `"`meth'"'=="irls_0" {
		local s1 "3" 
		if `"`2'"'!="" {
			capture confirm numeric variable `2'
			if _rc {
				di as err "init(`init0') invalid."
				confirm numeric variable `2'
				exit 198
			}
			else local mu0_var `"`2'"'
		}
		else {
			local mu0_var ""
		}
	}
	else {
		di as err "init(`init0') invalid. The available options for initial values are: irls_u, irls_r, irls_0 [varname]"
		exit 198
	}
	ret local start `s1'
	ret local mu0_var `mu0_var'
end
//-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// GLMWALS_MapFam: a modified version of MapFam in glm.ado, which maps the family
//-----------------------------------------------------------------------------------------------------------------------	
program GLMWALS_MapFam, rclass 
	version 14.0
	args f 
	local s1
	local f = lower(trim(`"`f'"'))
	local l = length(`"`f'"')
	if `"`f'"'=="" { 
		local s1 "1" 
	}
	else if `"`f'"'==substr("gaussian",1,max(`l',3))  { 
		local s1 "1" 
	}
	else if `"`f'"'==substr("normal",1,`l')           { 
		local s1 "1" 
	}
	else if `"`f'"'==substr("binomial",1,`l')         { 
		local s1 "2" 
	}
	else if `"`f'"'==substr("bernoulli",1,`l')        { 
		local s1 "2" 
	}
	else if `"`f'"'==substr("poisson",1,`l')          { 
		local s1 "3" 
	}
	else if `"`f'"'==substr("gamma",1,max(`l',3))     { 
		local s1 "4" 
	}
	else if `"`f'"'==substr("igaussian",1,max(`l',2)) { 
		local s1 "5" 
	}
	else if `"`f'"'==substr("inormal",1,max(`l',2))   { 
		local s1 "5" 
	}
	else if `"`f'"'=="ivg"                            { 
		local s1 "5" 
	}
	if "`s1'" != "" {
		local s1 "glmwals_v`s1'"
	}
	else {
		local s1 "`f'"
	}
	ret local famcode `s1'
end
//-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// GLMWALS_MapLink: a modified version of MapLink in glm.ado, which maps the link
//-----------------------------------------------------------------------------------------------------------------------	
program GLMWALS_MapLink, rclass 
	version 14.0 
	args ulink upow rest
	if `"`rest'"'!="" {
			di as err "link(`ulink' `upow' `rest') invalid"
			exit 198
	}
	local ulink = lower(`"`ulink'"')
	local l     = length(`"`ulink'"')
	local s1      "11"                /* power(p) link */
	local s2      .                   /* p of power(p) */
	if `"`ulink'"'=="" { 
		local s1 "" 
	}
	else if `"`ulink'"'==substr("identity",1,`l')   	{ 
		local s1 "01"
		local s2 1    
	}
	else if `"`ulink'"'==substr("reciprocal",1,`l') 	{ 
		local s1 "09"
		local s2 -1   
	}
    else if `"`ulink'"'=="log"                       	{ 
		local s1 "03"
		local s2 0    
	}
    else if `"`ulink'"'==substr("power",1,max(`l',3)) 	{
		capture confirm number `upow'
		if _rc {
				di as err "invalid # in link(power #)"
				exit 198
		}
		if `upow' == 0 {
			local s1 "03"
		}
		if `upow' == -1 {
			local s1 "09"
		}
		else if `upow' == -2 {
			local s1 "10"
		}
		local s2 `upow'
		local upow
	}
	else if `"`ulink'"'==substr("opower",1,max(`l',3)) {
		capture confirm number `upow'
		if _rc {
			di as err "invalid # in link(opower #)"
			exit 198
		}
		if `upow'==0 {
			local s1 "02"	
		}
		else {
			local s1 "12"   
		}
		local s2 `upow'
		local upow
	}
	else if `"`ulink'"'==substr("logit",1,`l')          { 
		local s1 "02"
		local s2 0    
	}
	else if `"`ulink'"'=="logc"    				        { 
		local s1 "05" 
	}
	else if `"`ulink'"'==substr("loglog",1,max(`l',4))  { 
		local s1 "06" 
	}
	else if `"`ulink'"'==substr("cloglog",1,`l')        { 
		local s1 "07" 
	}
	else if `"`ulink'"'==substr("probit",1,`l')         { 
		local s1 "08" 
	}
	else {
		local s1
		local s2 "`upow'"
	}
	if `"`s1'"' != "" {
		local s1 "glmwals_l`s1'"
	}
	else {
		local s1 "`ulink'"
	}
	ret local link  `s1'
	ret local power `s2'
end
//-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// GLMWALS_MapFL: a modified version of MapFL in glm.ado, which maps both the family and the link
//-----------------------------------------------------------------------------------------------------------------------	
program GLMWALS_MapFL, rclass 
	version 14.0
	args f ulink

	* Map the first argument of f
	GLMWALS_MapFam `f'					
	local fam `"`r(famcode)'"'			

	* Check and store (if any) the second argument of f (i.e. m for the bin family). 
	local mfixed 1
	local m 1
	if `"`fam'"'=="glmwals_v2" {           
		tokenize `"`f'"'
		if `"`2'"'!="" {
			capture confirm number `2'
			if _rc {
				unabbrev `2'
				local m `"`s(varlist)'"'
				local mfixed 0
			}
			else {
				if `2'>=.|`2'<1 {
					di as err `"`2' in family(binomial `2') invalid"'
					exit 198
				}
				local m `2'
			}
		}
		if `"`3'"'!="" {
			di as err "family(`f') invalid"
			exit 198
		}
	}
	
	* Store (if any) the second argument of a user-defined family. Additional arguments not allowed. 
	if `"`fam'"'~="glmwals_v1" & `"`fam'"'~="glmwals_v2" & /*
	*/ `"`fam'"'~="glmwals_v3" & `"`fam'"'~="glmwals_v4" & /*
	*/ `"`fam'"'~="glmwals_v5" {
		tokenize `"`f'"'
		if `"`2'"' != "" {
			global GLMWALS_fa `2'												
		}
		if `"`3'"'!="" {
			di as err "family(`f') invalid"
			exit 198
		}
	}

	* Map the link
	GLMWALS_MapLink `ulink'
	local link `"`r(link)'"'
	local pow  `"`r(power)'"'

	* If link is empty, then apply as default the canonical link of each family
	if `"`link'"'=="" {             
		local link "glmwals_l11"
		if `"`fam'"'=="glmwals_v1"      { 
			local link "glmwals_l01"
			local pow 1           
		}
		else if `"`fam'"'=="glmwals_v2" { 
			local link "glmwals_l02"
			local pow 0           
		}
		else if `"`fam'"'=="glmwals_v3" { 
			local link "glmwals_l03"
			local pow  0          
		}
		else if `"`fam'"'=="glmwals_v4" { 
			local link "glmwals_l09"
			local pow -1          
		}
		else if `"`fam'"'=="glmwals_v5" { 
			local link "glmwals_l10"
			local pow -2          
		}
		else local link
	}

	* Return locals
	ret local family `"`fam'"'
	ret local link   `"`link'"'
	ret local power  `pow'
	ret local m      `m'
	ret local mfixed `mfixed'
end
//-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// Di_glmwals: display of glmwals output 
//-----------------------------------------------------------------------------------------------------------------------	
program Di_glmwals
	syntax [, LEVel(cilevel) noTABle noFOCus noAUXiliary noHEADer  cformat(string) eform]
	local y=abbrev(`"`e(depvar)'"',10)
	local focvars 	"`e(focvars)'"
	local auxvars	"`e(auxvars)'"
	local k1=e(k1)
	local k2=e(k2)
	local regressors: coln e(b)
	local MLU_k		: colsof e(b)
	local offset `e(offset)'
	if "`offset'"!="" {
		local offrep: subinstr local offset "ln(" "", count(local offnum)
		if `offnum'>0 	local offtit "(exposure)"
		else 			local offtit "(offset)"
	}
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
		if `r(omit)'==0	 {
			local noncoll_X "`noncoll_X' `jj'"
		}
	}
	local LL=10
	foreach xx of local noncoll_X {
	    local ll: strlen local xx
		local LL=max(`ll', `LL')
	}
	if "`offset'"!="" {
		local ll: strlen local offset
		local LL=max(`ll', `LL')
	}
	tempname CI WALS_CI
	walsci `"`e(bcsimdata)'"' `level' `CI'
	matrix `WALS_CI'=J(2,`MLU_k',0)
	local ss=1
	forvalues jj=1(1)`MLU_k'{
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
				as text _col(`p1') "Number of obs =" 	
				as res %8.0f e(N)
			
			_n as text "Prior             : " as res "`e(prior)'"
			   as text _col(`p1') "Residual df " 	_col(`p2') "=" 	
			   as res %8.0f e(df_r) 							
			
			_n as text "Family            : " as res "`e(varfunct)'"
			   as text _col(`p1') "k1 " 			_col(`p2') "=" 	
			   as res %8.0f e(k1) 							
			
			_n as text "Link              : " as res "`e(linkt)'"
			   as text _col(`p1') "k2 " 			_col(`p2') "=" 	
			   as res %8.0f e(k2) 								
			
			_n as text "Variance function : " as res "V(u) = `e(varfuncf)'"
			   as text _col(`p1') "sigma" _col(`p2') "=" 
			   as res %8.0g e(sigma)
			
			_n as text "Link function     : " as res "g(u) = `e(linkf)'"
			   as text _col(`p1') "MC reps " 		_col(`p2') "=" 
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
		if "`eform'"=="" {
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
		}
		else {
			#delimit;
			di as text "{hline `LLp1'}{c TT}{hline 63}"													
				_n _col(`LL') "  {c |}"																
				_col(`p2') "`COEF_tit'"                            									
				_col(`p3') "`MOM_tit'"                            									
				_col(`p4') "`MOM_tit'"														 
				_col(`p5') "`MOM_tit'"														 
				_col(`p6') "`CI_tit'" 	                           									
				_n "{ralign `LL':`y'} {c |}"															
				"     exp(b)     Bias    Std.Err.    RMSE      [`level'% Conf. Int.]" 
				_n "{hline `LLp1'}{c +}{hline 63}";	
			#delimit cr	
		}
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
				*local xx_abbrev=abbrev("`xx'",15)
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
					if "`eform'"=="" {
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
					else {
						#delimit;
						di as text "{ralign `LL':`xx'} {c |}  "						
								_col(`p1') as res `cformat'  exp(_b[`xx']) 						"  "
								_col(`p2') as res `cformat'  `Bias'[1,`xx_num']*exp(_b[`xx'])	"  "
								_col(`p3') as res `cformat'  _se[`xx']*exp(_b[`xx'])			"  "
								_col(`p4') as res `cformat'  `RMSE'[1,`xx_num']*exp(_b[`xx'])	"  "
								_col(`p5') as res `cformat' exp(`CI'[1,`hh']) 					"  "
								_col(`p6') as res `cformat' exp(`CI'[2,`hh']);		
						#delimit cr	
					}
				}
				local hh=`hh'+1
			}
			local xx_num=`xx_num'+1
		}
		if "`offset'"!="" {
						#delimit;
						di as text "{ralign `LL':`offset'} {c |}  "						
								_col(`p1') as res `cformat'  1 				"  "
								_col(`p2') as text "`offtit'";		
						#delimit cr	
		}	
		noi di as text "{hline `LLp1'}{c BT}{hline 63}"		
	}
end
//-----------------------------------------------------------------------------------------------------------------------	

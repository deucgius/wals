* Weighted Average Least Squares (WALS) with heteroskedastic errors
*! Version 1.0
* Date			: 2024/10/05
* Authors		: De Luca Giuseppe, Jan R. Magnus
*-----------------------------------------------------------------------------------------------------------------------	
* hetwals
*-----------------------------------------------------------------------------------------------------------------------	
program define hetwals, eclass sort
	version 14.0, missing
	if replay() {
		if "`e(cmd)'" != "hetwals" {
			error 301
		}
		syntax [, SAVing(string asis) noHEADer noTABle noFOCus noAUXiliary LEVel(cilevel) cformat(string)]
		local bcsimdata `"`e(bcsimdata)'"'
		cap confirm file `"`bcsimdata'"'
		if _rc {
			di as err "MC replicates of bias-corrected hetwals estimates not found. You need to refit the model and resave the results" 
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
		Di_hetwals, `header' `table' `focus' `auxiliary' lev(`level') `cformat'
		exit
	}
	if !replay() {
		syntax [anything] [fw aw iw /] [if] [in], [*]
		local cmdline : copy local 0
	}
	Estimate `0'
	ereturn local cmdline `"hetwals `0'"'
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
syntax [anything] [fw aw iw /] [if] [in], 											
	het(varlist numeric min=1 ts fv)
	[														
	AUXCONstant												
	noCONstant												
	/* sigma(real 0)  not allowed */											
	
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
	
	TWOstep 
	HCONSTraints(string)
	waldhet 
	
	DIFficult 
	TECHnique(string) 
	ITERate(integer 300) 
	TOLerance(real 1e-6)
	LTOLerance(real 1e-7)
	NRTOLerance(real 1e-5) 
	NONRTOLerance
	from(string)
	SHOWHetreg
	noLOG

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
if "`:list dups het'"!="" {  
	di as err "{bf:het} cannot contain duplicate variables"
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
if `:list depvar in het' {
	di as err "dependent variable {bf:`depvar'} may not be included in {bf:het}"
	exit 198
}
if subinstr("`foclist0' `auxlist0'",".`depvar' ","",.) != "`foclist0' `auxlist0'" {
	di as err "time-series operators of dependent variable {bf:`depvar'} may not be " /*
	*/ 	"included as independent variables"
	exit 198
}
if subinstr("`het'",".`depvar' ","",.) != "`het'" {
	di as err "time-series operators of dependent variable {bf:`depvar'} may not be " /*
	*/ 	"included in {bf:het}"
	exit 198
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Mark estimation sample 
*------------------------------------------------------------------------------------------------
marksample touse
markout `touse' `depvar' `foclist0' `auxlist0' `exp' `het'
count if `touse' 
if r(N)==0 	error 2000 
if r(N)==1 	error 2001
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Options for heteroskedasticity
*------------------------------------------------------------------------------------------------
if "`twostep'" == "" {

	local mle "mle"
	
	* Options of hetregress, mle 
	#delimit;
	local hetopts "
		mle constraints(`hconstraints') `waldhet' `difficult' technique(`technique') iterate(`iterate')
		tolerance(`tolerance') ltolerance(`ltolerance') 
		nrtolerance(`nrtolerance') `nonrtolerance' from(`from') `log'
		"; 
	#delimit cr
	
	* Weights (allowed only with mle)
	if `"`weight'"' != "" {
		tempvar wgt_var
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
}	
else {	
	if "`hconstraints'" != "" {
		di as err "{bf:hconstraints} not allowed with {bf:twostep}"
		exit 198
	}
	if "`difficult'" != "" {
		di as err "{bf:difficult} not allowed with {bf:twostep}"
		exit 198
	}
	if "`technique'" != "" {
		di as err "{bf:technique} not allowed with {bf:twostep}"
		exit 198
	}
	if `iterate' != 300 {
		di as err "{bf:iterate} not allowed with {bf:twostep}"
		exit 198
	}
	if `tolerance' != 1e-6 {
		di as err "{bf:tolerance} not allowed with {bf:twostep}"
		exit 198
	}
	if `ltolerance' != 1e-7 {
		di as err "{bf:ltolerance} not allowed with {bf:twostep}"
		exit 198
	}
	if `nrtolerance' != 1e-5 {
		di as err "{bf:nrtolerance} not allowed with {bf:twostep}"
		exit 198
	}
	if "`nonrtolerance'" != "" {
		di as err "{bf:nonrtolerance} not allowed with {bf:twostep}"
		exit 198
	}
	if "`from'" != "" {
		di as err "{bf:from} not allowed with {bf:twostep}"
		exit 198
	}
	if "`log'"=="nolog" {
		di as err "{bf:nolog} not allowed with {bf:twostep}"
		exit 198
	}
	if `"`weight'"' != "" {
		noi di as error "weights are not allowed with option {bf:twostep}"
		exit 198
	}
	local hetopts "twostep"
}

* Saving
cap erase `"`c(tmpdir)'/STWALS_000001.tmp"'
noi _savingopt_parse CIMC_save CIMC_savereplace : saving ".dta" `"`saving'"'
if `"`CIMC_save'"'!="" & `"`CIMC_savereplace'"'=="" confirm new file `"`CIMC_save'"'
if `"`CIMC_save'"'=="" {
	local CIMC_save `"`c(tmpdir)'/STWALS_000001.tmp"'
	local CIMC_savereplace "replace"
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
if "`nc'"=="noconstant" & "`ac'"!="" {
	noi di in red "cannot specify both {bf:noconstant} and {bf:auxconstant}"
	error 184
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* GLS transformation for heteroskedastic
*------------------------------------------------------------------------------------------------
if "`showhetreg'"!="" {
    noi di as text "Estimates of unrestricted model"
	cap noi hetregress `depvar' `foclist0' `auxlist0' `wgt_exp' if `touse', `nc' het(`het') `hetopts'  
	if _rc!=0|("`e(method)'"=="ml" & "`e(converged)'"=="0") {
		noi "Convergence of heteroskedasticity model not achieved"
		error 430
	}	
}
else{
	noi _rmcoll `het' `wgt_exp' if `touse'
	cap hetregress `depvar' `foclist0' `auxlist0' `wgt_exp' if `touse', `nc' het(`het') `hetopts'  
	if _rc!=0|("`e(method)'"=="ml" & "`e(converged)'"=="0") {
		noi "Convergence of hetregress estimates for unrestricted model not achieved"
		error 430
	}	
}
replace `touse'=e(sample)
local N=e(N)
local het_method `e(method)'
local het_test=`e(chi2_c)'
local het_test_df=`e(df_m_c)'
local het_test_p=`e(p_c)'
local het_test_type "`e(chi2_ct)'"
tempname het_b het_V 
matrix `het_b'=e(b)
matrix `het_b'=`het_b'[1,"lnsigma2:"]
matrix `het_V'=e(V)
matrix `het_V'=`het_V'["lnsigma2:","lnsigma2:"]
tempvar het_sigma wgt_var_W
predict double `het_sigma' 					if `touse', sigma
gen double `wgt_var_W' = 1/`het_sigma'^2 	if `touse'
if "`twostep'"!="" {
	sum `wgt_var_W', meanonly
	replace `wgt_var_W'=`wgt_var_W'/r(mean)
	local wgt_exp_W `"[iw=`wgt_var_W']"'
	local sigma_W=0
}
else {
	sum `wgt_var_W' `wgt_exp', meanonly
	if `"`weight'"' != "" replace `wgt_var_W'=`wgt_var'*`wgt_var_W'/r(mean) if `touse'
	else 				  replace `wgt_var_W'=`wgt_var_W'/r(mean)			if `touse'
	local wgt_exp_W `"[iw=`wgt_var_W']"'
	regress `depvar' `foclist0' `auxlist0' `wgt_exp_W' if `touse'==1, `nc'
	local sigma_W=sqrt(e(rss)/e(N))
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* WALS estimates with rescaled iweights
*------------------------------------------------------------------------------------------------
if "`foclist0'"!="" {
	#delimit;
	noi wals `depvar' (`foclist0')  `auxlist0'  `wgt_exp_W' if `touse', 
		`nc' `ac' sigma(`sigma_W') prior(`PRI') `QUADRATURE' 
		plugin(`plugin') patht(`pathtab')
		lev(`level') rep(`reps') rseed(`rseed') saving(`saving') 
		notable nohead `fast';
	#delimit cr
}
else {
	#delimit;
	noi wals `depvar'  `auxlist0'  `wgt_exp_W' if `touse', 
		`nc' `ac' sigma(`sigma_W') prior(`PRI') `QUADRATURE' 
		plugin(`plugin') patht(`pathtab')
		lev(`level') rep(`reps') rseed(`rseed') saving(`saving') 
		notable nohead `fast';
	#delimit cr
}	
*------------------------------------------------------------------------------------------------




*------------------------------------------------------------------------------------------------
* Parse wals results
*------------------------------------------------------------------------------------------------
replace `touse'=e(sample)
local k0=e(k0) 
local k1=e(k1) 
local k2=e(k2)
local df_r=e(df_r)
local reps=e(reps)
local auxvars `e(auxvars)'
local focvars `e(focvars)'
local omitvars `e(omitvars)'
local allvars `e(allvars)'
local plugin `e(plugin)'
tempname WALS_b prior_par M_W M_pm M_t_ratios
matrix `WALS_b'=e(b)
matrix `prior_par'=e(priorpar)
matrix `M_W'=e(wals_wgt)
matrix `M_pm'=e(wals_pm)
matrix `M_t_ratios'=e(t_ratios)
if "`fast'"=="" {
	local k =`k1'+`k2'
	local wals_bc_fname `"`e(bcsimdata)'"'
	tempname WALS_bias WALS_V WALS_rmse WALS_MSE WALS_CI   
	matrix `WALS_bias'	=e(bias)
	matrix `WALS_V'		=e(V)
	matrix `WALS_rmse'	=e(rmse)
	matrix `WALS_MSE'	=e(MSE)
	matrix `WALS_CI'	=e(ci)
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Return estimation results 
*------------------------------------------------------------------------------------------------
if "`fast'"!="" {
	ereturn post `WALS_b' `wgt_exp_post', dep(`depvar') obs(`N') esample(`touse') dof(`df_r') buildfvinfo 
	ereturn scalar k0 			=`k0'
	ereturn scalar k1 			=`k1'
	ereturn scalar k2 			=`k2'
	ereturn scalar df_r			=`e(df_r)'				
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
	ereturn scalar hchi2		=`het_test'	
	ereturn scalar hchi2_df		=`het_test_df'	
	ereturn scalar hchi2_p		=`het_test_p'	

	
	ereturn local predict 		"hetwals_p"
	ereturn local properties	"`e(properties)'"
	if "`PRI'"=="laplace" 		ereturn local quadmethod	"analytic"
	else 						ereturn local quadmethod	"`quadm'"
	ereturn local prior 		"`PRI'"
	ereturn local wexp			"`e(wexp)'"
	ereturn local wtype			"`e(wtype)'"
	ereturn local htest			"`het_test_type'"	
	ereturn local hmethod		"`het_method'"	
	ereturn local hvars 		"`het'"
	if "`nc'"=="noconstant"						ereturn local constype "`nc'"
	else if "`nc'"~="noconstant" & "`ac'"=="" 	ereturn local constype "focus"
	else 										ereturn local constype "auxiliary"
	ereturn local auxvars 		"`auxvars'"
	ereturn local focvars 		"`focvars'"
	ereturn local omitvars	 	"`omitvars'"
	ereturn local allvars 		"`allvars'"
	ereturn local depvar		"`e(depvar)'"
	ereturn local title 		"Heteroskedastic WALS estimates"
	ereturn local cmd 			"hetwals"
	
	ereturn matrix priorpar		=`prior_par'
	ereturn matrix wals_wgt		=`M_W'
	ereturn matrix wals_pm		=`M_pm'
	ereturn matrix t_ratios		=`M_t_ratios'
	ereturn matrix het_V		=`het_V'	
	ereturn matrix het_b		=`het_b' 
}
else{
	ereturn post `WALS_b' `WALS_V' `wgt_exp_post', dep(`depvar') obs(`N') esample(`touse') dof(`df_r') buildfvinfo
 	ereturn scalar k0 			=`k0'
 	ereturn scalar k1 			=`k1'
 	ereturn scalar k2 			=`k2'
 	ereturn scalar rank 		=`k'
	ereturn scalar df_r			=`e(df_r)'				
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
 	ereturn scalar hchi2		=`het_test'	
 	ereturn scalar hchi2_df		=`het_test_df'	
 	ereturn scalar hchi2_p		=`het_test_p'	

 	ereturn local marginsok 	"XB default"
 	ereturn local predict 		"hetwals_p"
	ereturn local properties	"`e(properties)'"
 	ereturn local bcsimdata 	`"`wals_bc_fname'"'	
 	ereturn local plugin		"`plugin'"
 	if "`PRI'"=="laplace" 		ereturn local quadmethod	"analytic"
 	else 						ereturn local quadmethod	"`quadm'"
 	ereturn local prior 		"`PRI'"
	ereturn local wexp			"`e(wexp)'"
	ereturn local wtype			"`e(wtype)'"
 	ereturn local htest			"`het_test_type'"	
 	ereturn local hmethod		"`het_method'"	
 	ereturn local hvars 		"`het'"
 	if "`nc'"=="noconstant"						ereturn local constype "`nc'"
 	else if "`nc'"~="noconstant" & "`ac'"=="" 	ereturn local constype "focus"
 	else 										ereturn local constype "auxiliary"
 	ereturn local auxvars 		"`auxvars'"
 	ereturn local focvars 		"`focvars'"
 	ereturn local omitvars	 	"`omitvars'"
 	ereturn local allvars 		"`allvars'"
	ereturn local depvar		"`e(depvar)'"
 	ereturn local title 		"Heteroskedastic WALS estimates"
 	ereturn local cmd 			"hetwals"

 	ereturn matrix priorpar		=`prior_par'
 	ereturn matrix wals_wgt		=`M_W'
 	ereturn matrix wals_pm		=`M_pm'
 	ereturn matrix t_ratios		=`M_t_ratios'
 	ereturn matrix hV			=`het_V'	
 	ereturn matrix hb			=`het_b' 
 	ereturn matrix ci			=`WALS_CI'
 	ereturn matrix MSE			=`WALS_MSE'
 	ereturn matrix rmse			=`WALS_rmse'
 	ereturn matrix bias			=`WALS_bias'
}
*------------------------------------------------------------------------------------------------



*------------------------------------------------------------------------------------------------
* Display estimation results 
*------------------------------------------------------------------------------------------------
if "`fast'"=="" 	noi Di_hetwals, level(`level') `display_options'   
*------------------------------------------------------------------------------------------------
}
end
*-----------------------------------------------------------------------------------------------------------------------	



//-----------------------------------------------------------------------------------------------------------------------	
// Di_hetwals: display of wals output 
//-----------------------------------------------------------------------------------------------------------------------	
program Di_hetwals
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
	local chi : di %8.2f e(hchi2)
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
		#delimit;
		di as text _n `"`e(title)'"'					
			   as text _col(`p1') "Number of obs" 	_col(`p2') "="
			   as res %8.0f e(N)
			
			_n as text "Prior : " as res "`e(prior)'"
			   as text 	_col(`p1') "Residual df" 	_col(`p2') "=" 	
			   as res %8.0f e(df_r) 							
			
			_n as text "`e(htest)' test for heteroskedasticity"
			   as text _col(`p1') "k1 " 			_col(`p2') "=" 	
			   as res %8.0f e(k1) 							
			
			_n as text "chi2(" as res e(hchi2_df) as text ") = " as res `chi'
			   as text _col(`p1') "k2 " 			_col(`p2') "=" 	
			   as res %8.0f e(k2) 									
			
			_n as text "Prob > chi2 = " as res %6.4f e(hchi2_p)
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

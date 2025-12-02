*! Version 1.1
*Date: 2025/05/06
*Description: margins after wals estimation 
*Changes in the code: 
*    Version 1.1 - added message "Working on margins..." (line 173) 
*Author: De Luca Giuseppe, Jan R. Magnus
program define margwals, rclass
version 14.0, missing
qui {
	if 		"`e(cmd)'"!="wals"		&	/*
		*/	"`e(cmd)'"!="hetwals"	&	/*
		*/	"`e(cmd)'"!="ar1wals"	&	/*
		*/	"`e(cmd)'"!="xtwals"	&	/*
		*/	"`e(cmd)'"!="glmwals" {
		noi di as error  "wals was not the last estimation command"
		exit 301
	}

	* command line 
	local cmdline : copy local 0

	* Syntax 
	#delimit;	
	syntax [anything(name=marginlist)] [if] [in] [fw iw aw /], [
		
		PRedict(string asis)		
		dydx(string asis) 
		eyex(string asis)
		dyex(string asis)
		eydx(string asis)
		CONTinuous
		
		grand							
		
		at(string asis) 
		atmeans
		ASBALanced
		
		over(string asis) 				
		/* subpop(string asis) 			- EXCLUDED	*/
		within(string asis) 			
		/* contrast_options				- EXCLUDED	*/
		/* pwcompare_options			- EXCLUDED	*/
		
		/* vce(delta)					- ALWAYS	*/
		/* vce(unconditional) 			- EXCLUDED	*/
		/* nose							- EXCLUDED	*/ 
		
		NOWEIGHTs 					
		NOEsample
		emptycells(string asis)			 
		ESTIMTOLerance(string asis) 
		noestimcheck 
		force
		NOCHAINrule						
		CHAINrule						
		
		/* level(string asis) */
		LEVel(cilevel)
		/* mcompare(method) 			- EXCLUDED 	*/ 
		/* noatlegend					- ALWAYS 	*/	 
		/* post							- ALWAYS	*/
		/* display_options				- EXCLUDED	*/
		/* df(#) 						- EXCLUDED	*/
		];
	#delimit cr	

	* Margins options 
	if "`predict'"!=""			local dydx "predict(`predict')"
	if "`dydx'"!=""				local dydx "dydx(`dydx')"
	if "`eyex'"!=""				local eyex "eyex(`eyex')"
	if "`dyex'"!=""				local dyex "dyex(`dyex')"
	if "`eydx'"!=""				local eydx "eydx(`eydx')"
	if "`at'"!=""				local at "at(`at')"
	if "`over'"!=""				local over "over(`over')"
	if "`within'"!=""			local within "within(`within')"
	if "`emptycells'"!=""		local emptycells "emptycells(`emptycells')"
	if "`estimtolerance'"!=""	local estimtolerance "estimtolerance(`estimtolerance')"
	#delimit;
		local mar_opts 
				"`predict' `dydx' `eyex' `dyex' `eydx' 
				`continuous' `grand' `at' `atmeans' `asbalanced'
				`over' `within'
				`noweights' `noesample'  `emptycells' `estimtolerance' 
				`force' `chainrule' `nochainrule'  
				noestimcheck noatlegend nopvalue nose";
	#delimit cr
	local mar_opts: list retokenize mar_opts

	* Weights
	if `"`weight'"' != "" {
		if `"`weight'"'=="pweight" {
			noi di as err "`weight' not allowed"
			error 198
		}
		else if `"`weight'"'=="aweight" {
			local wgt_exp `"[aw=`exp']"'		
		}
		else if `"`weight'"'=="fweight" {
			local wgt_exp `"[fw=`exp']"'
		}
		else {
			local wgt_exp `"[iw=`exp']"'
		}
	}
 
	* Latest wals estimates and associated objects   
	tempname wals_est 
	estimates store `wals_est'

	* Point estimates of WALS margins 
	margins `marginlist' `if' `in' `wgt_exp', `mar_opts'
	tempname wals_mar_est wals_mar_b
	_return hold `wals_mar_est'
	_return restore `wals_mar_est', hold
	matrix `wals_mar_b'=r(b)
	local wals_mar_n0: colsof `wals_mar_b'
	local wals_mar_names0: colfullnames `wals_mar_b'

	* list and number of all covariates 
	local Xlist0 "`e(allvars)'"	
	local k0: word count `e(allvars)'

	* list and number of non-collinear covariates 
	local Xlist "`e(focvars)' `e(auxvars)'"				
	if "`e(constype)'"!="noconstant" {
		local Xlist: subinstr local Xlist "_cons" ""	
		local Xlist "`Xlist' _cons"	
	}
	local k=e(rank)

	* list and number of omitted covariates 
	local Xlist_omit: list Xlist0 - Xlist

	* list of variables in the simulation dataset
	local wals_bc_list ""
	forvalues h=1(1)`k'{
		local wals_bc_list "`wals_bc_list' wals_bc_`h'"
	}

	* Parse to MATA extended dataset of MC replications (including omitted covariates)
	preserve 
	cap use `"`e(bcsimdata)'"', clear
	if _rc!=0 {
		noi di as err "{p 4 4 2}MC replications of bias-corrected wals estimates, "	/*
			*/ "as previously saved by {bf:`e(cmd)'}, do not exist. You need "	/*
			*/ "to refit the model.{p_end}" 
		exit _rc	
	}
	local rwals_bc_list ""
	local oo=1
	local ss=1
	forvalues rr=1(1)`k0'{
		local xx: word `rr' of `Xlist0'
		local xo: word `oo' of `Xlist_omit'
		if "`xx'"!="`xo'" {
			rename wals_bc_`ss' rwals_bc_`rr'
			local ++ss
		}
		else {
			gen double rwals_bc_`rr'=0
			lab var rwals_bc_`rr' "`xo'"
			local ++oo
		}
		local rwals_bc_list "`rwals_bc_list' rwals_bc_`rr'"
	}
	order `rwals_bc_list'
	tempname wals_bc_repdata
	mata: `wals_bc_repdata'=st_data(., "`rwals_bc_list'")
	restore

	* MC replications of bias-corrected WALS margins
	noi di in gr "Working on margins..." 
	local wals_reps `e(reps)'
	noi margwals_reps `marginlist' `if' `in' `wgt_exp', /*
		*/ bc_repdata(`wals_bc_repdata') 				/*
		*/ mar_n0(`wals_mar_n0') mar_opts(`mar_opts')	/*
		*/ reps(`wals_reps') 
	
	local wals_mar_repdata `e(mar_repdata)'

	* Store MC replications of bias-corrected WALS margins
	local MAR_save `"`c(tmpdir)'/STWALS_000002.tmp"'
	preserve 
		clear
		set obs `wals_reps'
		forvalues h=1(1)`wals_mar_n0' {
		    local ll: word `h' of `wals_mar_names0'
			gen double wals_bc_`h'=.
			label var wals_bc_`h' "`ll'"
		}
		mata: st_store(.,., `wals_mar_repdata')
		local wals_mar_names
		local wals_mar_n=1
		forvalues h=1(1)`wals_mar_n0' {
		    local ll: word `h' of `wals_mar_names0'
			_ms_parse_parts `ll'
			if `r(omit)'==1 {
			    drop wals_bc_`h'
			}
			else {
			    rename wals_bc_`h' wals_bc_`wals_mar_n'
				local wals_mar_names "`wals_mar_names' `ll'"
			    local wals_mar_n=`wals_mar_n'+1
			}
		}
		local wals_mar_names=trim(`"`wals_mar_names'"')
		local wals_mar_n=`wals_mar_n'-1
		save `"`MAR_save'"', replace
	restore

	* Confidence intervals of wals margins
	tempname CI wals_mar_CI 
	walsci `"`MAR_save'"' `level' `CI'
	matrix `wals_mar_CI'=J(2,`wals_mar_n0',0)
	local ss=1
	forvalues jj=1(1)`wals_mar_n0'{
		local xx: word `jj' of `wals_mar_names0'
		local yy: word `ss' of `wals_mar_names'
		if "`xx'"=="`yy'" {
			matrix `wals_mar_CI'[1,`jj']=`CI'[1,`ss']
			matrix `wals_mar_CI'[2,`jj']=`CI'[2,`ss']
			local ++ss
		}
	}
	matrix rown `wals_mar_CI'= ll ul
	matrix coln `wals_mar_CI'=`wals_mar_names0'
	cap erase `"`c(tmpdir)'/STWALS_000002.tmp"'	
	
	* Restore wals estimates
	estimates restore `wals_est'

	* Other results from margins  
	_return restore `wals_mar_est', hold

	local N 			= r(N)
	local k_margins		= r(k_margins)
	local k_predict		= r(k_predict)
	local k_at			= r(k_at)
	
	local title 		"`r(title)'"
	local expression	"`r(expression)'"
	local est_cmd		"`r(est_cmd)'"
	local est_cmdline	"`r(est_cmdline)'"
	local xvars			"`r(xvars)'"
	local margins		"`r(margins)'"
	forvalues jj=1(1)`k_predict' {
		local predict`jj'_opts 	"`r(predict`jj'_opts)'"
		local predict`jj'_label "`r(predict`jj'_label)'"
	}
	local derivatives 	"`r(derivatives)'"
	local emptycells	"`r(emptycells)'"
	local over			"`r(over)'"
	local within		"`r(within)'"
	
	tempname uN Err
	matrix `uN'			= r(_N)
	matrix `Err'		= r(error)
	if `k_at'>0 {
	    tempname at
		matrix `at'		= r(at)
		local n_at: rowsof `at'
		forvalues jj=1(1)`n_at' {
		    local atstats`jj'	"`r(atstats`jj')'"
		}
	}

	* Equation type 
	local eq_list: coleq `wals_mar_b'
	local eq_1: word 1 of `eq_list'
	if "`eq_1'"=="_" 	local eq_type "empty"
	else 				local eq_type "full"	

	* Coef title
	if 	 "`title'"=="Predictive margins"			///
		|"`title'"=="Adjusted predictions" 	{
		    local COEF="Margin"
	}
	else {
	    local COEF=" `derivatives'"
	}

	* Display header
	noi di as text "`title'"  _col(40) "Number of obs = " as res `N' /*
			*/ _n as text 	  _col(40) "MC reps       = " as res `wals_reps'
			
	* Display coefficient table
	local LL=22
	foreach xx of local wals_mar_names {    
		local ll: strlen local xx
		local LL=max(`LL', `ll')
	}
	foreach xx of local eq_list {
		local ll: strlen local xx
		local LL=max(`LL', `ll')
	}
	if "`e(plugin)'"=="ds" {
		local CI_tit 	"MC-DS"
	}
	else {
		local CI_tit 	"MC-ML"
	}
	local LLp1=`LL'+1
	local p2=`LL'+8	
	local p3=`LL'+27
	#delimit;
	noi di as text "{hline `LLp1'}{c TT}{hline 36}"													
		_n _col(`LL') "  {c |}"																
		_col(`p2') "`COEF_tit'"                            									
		_col(`p3') "`CI_tit'" 	                           									
		_n "{ralign `LL':`y'} {c |}"															
		"     Coef.      [`level'% Conf. Interval]" 
		_n "{hline `LLp1'}{c +}{hline 36}";	
	#delimit cr	

	if "`eq_type'"=="empty" {
		local p1=`LL'+3	
		local p2=`LL'+18	
		local p3=`LL'+30
		local hh=1	
		foreach xx of local wals_mar_names0 {
			_ms_parse_parts `xx'
			if `r(omit)'==0 & (`wals_mar_b'[1,`hh']!=0|`wals_mar_CI'[1,`hh']!=0|`wals_mar_CI'[2,`hh']){
				#delimit;
				noi di as text "{ralign `LL':`xx'} {c |}  "						
						_col(`p1') as res %9.0g  `wals_mar_b'[1,`hh'] 			"  "
						_col(`p2') as res %9.0g  `wals_mar_CI'[1,`hh']			"  "
						_col(`p3') as res %9.0g  `wals_mar_CI'[2,`hh'];
				#delimit cr	
			}
			local hh=`hh'+1
		}
		noi di as text "{hline `LLp1'}{c BT}{hline 36}"		
	}
	else {
		local p1=`LL'+3	
		local p2=`LL'+18	
		local p3=`LL'+30
		local hh=1	
		local eq_name_last ""
		foreach xx of local wals_mar_names0 {
			_ms_parse_parts `xx'
			if `r(omit)'==0 & (`wals_mar_b'[1,`hh']!=0|`wals_mar_CI'[1,`hh']!=0|`wals_mar_CI'[2,`hh']){

				* Equation name
				local eq_name: word `hh' of `eq_list'
				if "`eq_name'"!="`eq_name_last'" {
					if "`eq_name_last'"!="" noi di as text "{hline `LLp1'}{c +}{hline 36}"
					noi di as text "`eq_name'" _col(`LL') "  {c |}"
					local eq_name_last "`eq_name'"
				}
				
				* Estimated effects of levels within secondary names
				#delimit;
				noi di as text "{ralign `LL':`xx'} {c |}  "						
						_col(`p1') as res %9.0g `wals_mar_b'[1,`hh'] 			"  "
						_col(`p2') as res %9.0g `wals_mar_CI'[1,`hh'] 			"  "
						_col(`p3') as res %9.0g `wals_mar_CI'[2,`hh'];		
				#delimit cr	
			}
			local hh=`hh'+1
		}
		noi di as text "{hline `LLp1'}{c BT}{hline 36}"		
	}
	
	* Return results
	return scalar level		=`level'
	return scalar k_at		=`k_at'
	return scalar k_predict =`k_predict'
	return scalar k_margins =`k_margins'
	return scalar N			=`N'
	
	return matrix ci		=`wals_mar_CI'
	return matrix b			=`wals_mar_b'
	if `k_at'>0 {
	    return matrix at 	=`at'
	}
	return matrix error		=`Err'
	return matrix _N		=`uN'
	
	return local emptycells		"`emptycells'"
	return local within			"`within'"
	return local over			"`over'"
	if `k_at'>0 {
		forvalues jj=1(1)`n_at' {
			return local atstats`jj' "`atstats`jj''"
		}
	}
	return local derivatives	"`derivatives'"
	return local xvars			"`xvars'"
	forvalues jj=1(1)`k_predict' {
		return local predict`jj'_opts  "`predict`jj'_opts'"
		return local predict`jj'_label "`predict`jj'_label'"
	}
	return local expression		"`expression'"
	return local title			"`title'"
	return local est_cmd		"`est_cmd'"
	return local est_cmdline	"`est_cmdline'"
	return local cmd 			"margwals" 
	return local cmdline 		"margwals`cmdline'"	
}					
end
//-----------------------------------------------------------------------------------------------------------------------	


//-----------------------------------------------------------------------------------------------------------------------	
// margwals_reps
//-----------------------------------------------------------------------------------------------------------------------	
program define margwals_reps, eclass
version 14.0, missing
qui {
	syntax [anything(name=marginlist)] [if] [in] [fw iw aw /], 	/*
		*/ bc_repdata(string asis) mar_n0(integer) 				/*
		*/ mar_opts(string asis) reps(integer) 	

	if `"`weight'"' != "" {
		if `"`weight'"'=="pweight" {
			noi di as err "`weight' not allowed"
			error 198
		}
		else if `"`weight'"'=="aweight" {
			local wgt_exp `"[aw=`exp']"'		
		}
		else if `"`weight'"'=="fweight" {
			local wgt_exp `"[fw=`exp']"'
		}
		else {
			local wgt_exp `"[iw=`exp']"'
		}
	}
	tempname mar_repdata wals_bc_jj
	mata: `mar_repdata'=J(0,`mar_n0',0)
	forvalues jj=1(1)`reps' {
		mata: st_matrix("`wals_bc_jj'", `bc_repdata'[`jj',.])
		ereturn repost b = `wals_bc_jj' , properties(b) buildfvinfo 
		margins `marginlist' `if' `in' `wgt_exp', `mar_opts'
		mata: `mar_repdata'=(`mar_repdata' \ st_matrix("r(b)"))
	}	
	ereturn local mar_repdata "`mar_repdata'"
}
end
//-----------------------------------------------------------------------------------------------------------------------	



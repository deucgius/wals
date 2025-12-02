*! version 1.1
*! predictions of xtwals
*! De Luca Giuseppe, Jan R. Magnus (2024/10/05)
program define xtwals_p, sort
	if "`e(cmd)'" != "xtwals" {
		noi di as error  "xtwals was not the last estimation command"
		exit 301
	}
	version 14.0
	#delimit; 
		syntax [anything] [if] [in] [,
			xb biasp stdp rmsep bcp stdbcp 
			ie bcie Residuals BCResiduals wp bcwp 
			CINTerval PINTerval LEVel(cilevel) 
			rseed(numlist max=1 >0 integer)];
	#delimit cr
	if "`cinterval'"=="" & "`pinterval'"=="" {	
		
		* Parse options of the first syntax (except standard options list xb)
		local myopts "biasp stdp rmsep bcp stdbcp ie bcie Residuals BCResiduals wp bcwp"
		noi _pred_se "`myopts'" `0'
		if `s(done)' { 
			exit
		}
		local vtyp `s(typ)'
		local varn `s(varn)'
		local 0    `"`s(rest)'"'
		syntax [if] [in] [, `myopts']
		local type "`biasp'`stdp'`rmsep'`bcp'`stdbcp'`ie'`bcie'`residuals'`bcresiduals'`wp'`bcwp'"
		
		* default: xb assumed 
		if ("`type'"=="") {
			di in gr "(option xb assumed; fitted values)"
			tempvar lp
			qui _predict double `lp' `if' `in', xb
			gen `vtyp' `varn' = `lp' `if' `in'
			label var `varn' "Linear prediction"
			exit
		}
		* bias of linear prediction 
		else if "`type'"=="biasp" {
			qui {
				tempname bias 
				matrix `bias'=e(bias)
				local Xlist: coln e(b)
				local ctype "`e(constype)'"
				tempvar bp
				local k:colsof e(bias)
				gen double `bp'=`bias'[1,`k'] 	`if' `in'
				local h 1
				foreach v of local Xlist {
					if "`v'"!="_cons" qui replace `bp'=`bp'+`bias'[1,`h']*`v'
					local h=`h'+1
				}
			}
			gen `vtyp' `varn' = `bp' `if' `in'
			label var `varn' "Bias of the prediction"
			exit
		}
		* SE of linear prediction 
		else if "`type'"=="stdp" {
			tempvar sep
			qui _predict double `sep' `if' `in', stdp
			gen `vtyp' `varn' = `sep' `if' `in'
			label var `varn' "SE of the prediction"
			exit
		}
		* RMSE of linear prediction 
		else if "`type'"=="rmsep" {
			qui {
				tempvar bp sep
				predict double `bp'   `if' `in', biasp
				_predict double `sep' `if' `in', stdp
			}
			gen `vtyp' `varn' = sqrt(`bp'^2+`sep'^2) `if' `in'
			label var `varn' "RMSE of the prediction"
			exit
		}
		* bias-corrected linear prediction 
		else if "`type'"=="bcp" {
			qui {
				tempvar lp bp
				_predict double `lp' `if' `in', xb
				predict double `bp'   `if' `in', biasp
			}
			gen `vtyp' `varn' = `lp'-`bp' `if' `in'
			label var `varn' "Bias-corrected prediction"
			exit
		}
		* SE of bias-corrected linear prediction 
		else if "`type'"=="stdbcp" {
			qui{
				local Xlist "`e(focvars)' `e(auxvars)'"				
				local Xlist: subinstr local Xlist "_cons" ""	
				local Xlist "`Xlist' _cons"	
				preserve 
					tempname Vs
					use `"`e(bcsimdata)'"', clear
					corr, cov
					matrix `Vs'=r(C)
				restore
				tempvar const
				qui gen double `const'=1		`if' `in'
				local cons _cons
				local Xlist2: list Xlist - cons
				local Xlist2 "`Xlist2' `const'"
				fvrevar `Xlist2'
				local Xlist2 "`r(varlist)'"
				
				marksample touse
				gsort -`touse'
				count if `touse'==1
				local N=r(N)
				tempvar bcs
				gen double `bcs'=.	
				tempname xj xjVsxjt
				forvalues j=1(1)`N' {
					mkmat `Xlist2' if _n==`j', matrix(`xj')
					matrix `xjVsxjt'=`xj' * `Vs' * `xj''
					replace `bcs'=sqrt(`xjVsxjt'[1,1]) if _n==`j'
				}
			}
			gen `vtyp' `varn' = `bcs' if `touse'
			label var `varn' "S.E. of bias-corrected prediction"
			exit
		}
		* individual effects (in-sample only)
		else if "`type'"=="ie" {
			qui {
				* Fixed-effects
				if "`e(model)'"=="fe-iid" {
					tempvar esample ESAMPLE lp ie IE
					gen byte `esample'=e(sample)
					bys `e(ivar)': egen byte `ESAMPLE'=max(`esample')
					_predict double `lp' if `esample', xb
					capture xtset
					if "`r(timevar)'" != "" {
						tsrevar `e(depvar)'
						local depvar `r(varlist)'
					}
					else {
						local depvar `e(depvar)'
					}
					sort `e(ivar)' `esample'
					by `e(ivar)' `esample': gen double `ie' = sum(`depvar')/_n - sum(`lp')/_n 	if `esample'
					by `e(ivar)' `esample': replace `ie' = `ie'[_N] 							if `esample'
					by `e(ivar)': egen double `IE' = mean(`ie')  								if `ESAMPLE'
				}
				* Random-effects
				if "`e(model)'"=="re-iid" {
					tempvar esample ESAMPLE lp ie IE
					gen byte `esample'=e(sample)
					bys `e(ivar)': egen byte `ESAMPLE'=max(`esample')
					_predict double `lp' if `esample', xb

					capture xtset
					sort `e(ivar)' `esample'
					if "`r(timevar)'" != "" {
						tsrevar `e(depvar)'
						local depvar `r(varlist)'
					}
					else {
						local depvar `e(depvar)'
					}
					if e(Tcon)==1 {
						tempname ratio
						scalar `ratio' = scalar(e(sig_nu)^2/(e(Tavg)*e(sig_nu)^2+e(sig_e)^2))
					}
					else {
						tempvar ratio
						by `e(ivar)' `esample': gen double `ratio'=scalar(e(sig_nu)^2/(_N*e(sig_nu)^2+e(sig_e)^2)) if `esample'
					}
					sort `e(ivar)' `esample'
					by `e(ivar)': gen double `ie' = `ratio'*sum(`depvar'-`lp') 	if `esample'
					by `e(ivar)': replace `ie' = `ie'[_N] 						if `esample'
					by `e(ivar)': egen double `IE' = mean(`ie')  				if `ESAMPLE'
				}
				if "`e(errors)'"=="ar1" {
					noi di in red "option {bf:ie} not allowed in models with AR(1) errors"
					exit 198
				}
			}
			gen `vtyp' `varn' = `IE' 	`if' `in'
			label var `varn' "individual effects"
			exit
		}
		* bias-corrected individual effects (in-sample only)
		else if "`type'"=="bcie" {
			qui {
				* Fixed-effects iid
				if "`e(model)'"=="fe-iid" {
					tempvar esample ESAMPLE bcp bcie bcIE
					gen byte `esample' = e(sample)
					bys `e(ivar)': egen byte `ESAMPLE'=max(`esample')
					predict double `bcp' if `esample', bcp
					capture xtset
					if "`r(timevar)'" != "" {
						tsrevar `e(depvar)'
						local depvar `r(varlist)'
					}
					else {
						local depvar `e(depvar)'
					}
					sort `e(ivar)' `esample'
					by `e(ivar)' `esample': gen double `bcie' = sum(`depvar')/_n - sum(`bcp')/_n 	if `esample'
					by `e(ivar)' `esample': replace `bcie' = `bcie'[_N] 							if `esample'
					by `e(ivar)': egen double `bcIE' = mean(`bcie')  								if `ESAMPLE'
				}
				* Random-effects iid
				if "`e(model)'"=="re-iid" {
					tempvar esample ESAMPLE bcp bcie bcIE
					gen byte `esample' = e(sample)
					bys `e(ivar)': egen byte `ESAMPLE'=max(`esample')
					predict double `bcp' if `esample', bcp

					capture xtset
					sort `e(ivar)' `esample'
					if "`r(timevar)'" != "" {
						tsrevar `e(depvar)'
						local depvar `r(varlist)'
					}
					else {
						local depvar `e(depvar)'
					}
					sort `e(ivar)' `esample'
					if e(Tcon)==1 {
						tempname ratio
						scalar `ratio' = scalar(e(sig_nu)^2/(e(Tavg)*e(sig_nu)^2+e(sig_e)^2))
					}
					else {
						tempvar ratio
						by `e(ivar)' `esample': gen double `ratio'=scalar(e(sig_nu)^2/(_N*e(sig_nu)^2+e(sig_e)^2)) if `esample'
					}
					sort `e(ivar)' `esample'
					by `e(ivar)': gen double `bcie' = `ratio'*sum(`depvar'-`bcp') 	if `esample'
					by `e(ivar)': replace `bcie' = `bcie'[_N] 						if `esample'
					by `e(ivar)': egen double `bcIE' = mean(`bcie')  				if `ESAMPLE'					
				}
				if "`e(errors)'"=="ar1" {
					noi di in red "option {bf:bcie} not allowed in models with AR(1) errors"
					exit 198
				}
			}
			gen `vtyp' `varn' = `bcIE' 	`if' `in'
			label var `varn' "bias-corrected individual effects"
			exit
		}
		* residuals
		else if "`type'"=="residuals" {
			qui {
				if "`e(errors)'"=="iid" { 
					tempvar esample lp ie
					gen byte `esample' = e(sample)
					_predict double `lp' 	if `esample', xb
					predict double `ie' 	if `esample', ie
				}
				if "`e(errors)'"=="ar1" {
					noi di in red "option {bf:residuals} not allowed in models with AR(1) errors"
					exit 198
				}
			}
			gen `vtyp' `varn' = `e(depvar)' - `lp' - `ie' `if' `in'
			label var `varn' "residuals"
			exit			
		}
		* bias-corrected residuals
		else if "`type'"=="bcresiduals" {
			qui {
				if "`e(errors)'"=="iid" { 
					tempvar esample bcp bcie
					gen byte `esample' = e(sample)
					predict double `bcp' 	if `esample', bcp
					predict double `bcie' 	if `esample', bcie
				}
				if "`e(errors)'"=="ar1" {
					noi di in red "option {bf:bcresiduals} not allowed in models with AR(1) errors"
					exit 198
				}
			}
			gen `vtyp' `varn' = `e(depvar)' - `bcp' - `bcie' `if' `in'
			label var `varn' "bias-corrected residuals"
			exit			
		}
		* s periods ahead WALS predictions 
		else if "`type'"=="wp" {
			if "`e(errors)'"=="iid" { 
				qui {
					* Estimation sample
					tempvar esample ESAMPLE 
					gen byte `esample'=e(sample)
					bys `e(ivar)': egen byte `ESAMPLE'=max(`esample')
					
					* List of non-collinear regressors 
					local Xlist "`e(focvars)' `e(auxvars)'"				
					local Xlist: subinstr local Xlist "_cons" ""	
						
					* Mark sample (out-sample only)
					marksample touse
					markout `touse' `Xlist'
					replace `touse'=0 if `esample'==1
					replace `touse'=0 if `ESAMPLE'==0
					
					* Linear prediction
					tempvar lp
					_predict double `lp' if `touse'==1, xb
					
					* Individual effects
					tempvar IE 			
					predict double `IE' if `ESAMPLE'==1, ie
					replace `IE'=. 		if `touse'==0
				}
				if "`e(errors)'"=="ar1" {
					noi di in red "option {bf:wp} not allowed in models with AR(1) errors"
					exit 198
				}
			}
			gen `vtyp' `varn' = `lp'+`IE' 	if `touse'
			label var `varn' "s periods ahed WALS predictions"
			exit
		}
		* s periods ahead bias-corrected WALS predictions 
		else if "`type'"=="bcwp" {
			if "`e(errors)'"=="iid" { 
				qui {
					* Estimation sample
					tempvar esample ESAMPLE 
					gen byte `esample'=e(sample)
					bys `e(ivar)': egen byte `ESAMPLE'=max(`esample')
					
					* List of non-collinear regressors 
					local Xlist "`e(focvars)' `e(auxvars)'"				
					local Xlist: subinstr local Xlist "_cons" ""	
						
					* Mark sample (out-sample only)
					marksample touse
					markout `touse' `Xlist'
					replace `touse'=0 	if `esample'==1
					replace `touse'=0 	if `ESAMPLE'==0

					* Linear prediction
					tempvar bcp
					predict double `bcp' 	if `touse'==1, bcp
					
					* Individual effects
					tempvar bcIE 			
					predict double `bcIE' 	if `ESAMPLE'==1, bcie
					replace `bcIE'=. 	 	if `touse'==0
				}
				if "`e(errors)'"=="ar1" {
					noi di in red "option {bf:wp} not allowed in models with AR(1) errors"
					exit 198
				}
			}
			gen `vtyp' `varn' = `bcp'+`bcIE' 	if `touse'
			label var `varn' "s periods ahed bias-corrected WALS predictions"
			exit
		}
	
	}
	else {
		qui {
			* Parse options of the second syntax
			local myopts "CINTerval PINTerval LEVel(cilevel) rseed(numlist max=1 >0 integer)"
			gettoken bfc 0 : 0, parse(",")
			local ifin: list bfc - anything 
			local narg: word count `anything'
			if `narg'!=2 & `narg'!=3 {
				error 198
			}
			if `narg'==2 {
				local vtyp 
				local varn_lb: word 1 of `anything'
				local varn_ub: word 2 of `anything'
				local 0    "`ifin'`0'"
				syntax [if] [in] [, `myopts']
			}
			if `narg'==3 {
				local vtyp	 : word 1 of `anything'
				local varn_lb: word 2 of `anything'
				local varn_ub: word 3 of `anything'
				noi _pred_se "`myopts'" `vtyp' `varn_lb' `ifin'`0'
				local vtyp `s(typ)'
				local 0    `"`s(rest)'"'
				syntax [if] [in] [, `myopts']
			}
			cap opts_exclusive "`cinterval' `pinterval'"
			if _rc {
				noi di as err "only one of {bf:cinterval} or {bf:pinterval} is allowed"
				exit _rc
			}
			local type "`cinterval'`pinterval'"
			if "`e(errors)'"=="ar1" {
				noi di in red "option {bf:`type'} not allowed in models with AR(1) errors"
				exit 198
			}
			cap confirm new variable `varn_lb', exact 
			if _rc!=0 {
				noi di in red "variable `varn_lb' already defined"
				exit 110
			}
			cap confirm new variable `varn_ub', exact 
			if _rc!=0 {
				noi di in red "variable `varn_ub' already defined"
				exit 110
			}

			* Estimation sample
			tempvar esample ESAMPLE tvar
			gen byte `esample'=e(sample)
			bys `e(ivar)': egen byte `ESAMPLE'=max(`esample')
			bys `e(ivar)': gen long `tvar'=sum(`ESAMPLE') if `ESAMPLE'
			sort `e(ivar)' `tvar'

			* List of non-collinear covariates 
			local Xlist "`e(focvars)' `e(auxvars)'"				
			local Xlist: subinstr local Xlist "_cons" ""	
			tempname const
			gen double `const'=1 							
			local Xlist "`Xlist' `const'"	
							
			* Mark sample (out-sample only)
			marksample touse
			markout `touse' `Xlist'
			replace `touse'=0 if `esample'==1
			replace `touse'=0 if `ESAMPLE'==0
			count if `touse'==1
			local N=r(N)

			* List of regressors  (as temp vars)
			fvrevar `Xlist'							
			local M_vars "`r(varlist)'"

			* depvar
			capture xtset
			if "`r(timevar)'" != "" {
				tsrevar `e(depvar)'
				local depvar `r(varlist)'
			}
			else {
				local depvar `e(depvar)'
			}
			
			* Individual means
			tempvar ybar_i YBAR_i
			sort `e(ivar)' `esample' `tvar'
			by `e(ivar)' `esample': gen double `ybar_i' = sum(`depvar')/_n 	if `esample'
			by `e(ivar)': replace `ybar_i' = `ybar_i'[_N] 					if `esample'
			by `e(ivar)': egen double `YBAR_i' = mean(`ybar_i') 			if `ESAMPLE'
			replace `YBAR_i'=. 												if `touse'==0
			
			local h 1
			local xbar_list 
			local XBAR_list
			foreach v of local M_vars {
				tempvar xbar_`h' XBAR_`h' 
				sort `e(ivar)' `esample' `tvar'
				by `e(ivar)' `esample': gen double `xbar_`h'' = sum(`v')/_n if `esample'
				by `e(ivar)': replace `xbar_`h'' = `xbar_`h''[_N] 			if `esample'
				by `e(ivar)': egen double `XBAR_`h'' = mean(`xbar_`h'') 	if `ESAMPLE'
				replace `XBAR_`h''=. 										if `touse'==0
				local xbar_list "`xbar_list' `xbar_`h''"
				local XBAR_list "`XBAR_list' `XBAR_`h''"
				local h=`h'+1
			}

			* Estimates of variance components (only re model)
			if e(Tcon)==1 {
				tempname ratio
				scalar `ratio' = scalar(e(sig_nu)^2/(e(sig_nu)^2+e(sig_e)^2/e(Tavg)))
			}
			else {
				tempvar ratio RATIO
				by `e(ivar)' `esample': gen double `ratio'=scalar(e(sig_nu)^2/(e(sig_nu)^2+e(sig_e)^2/_N)) 	if `esample'
				by `e(ivar)': egen double `RATIO' = mean(`ratio') 											if `ESAMPLE'
				replace `RATIO'=. 																			if `touse'==0
			}
					
			* MC replications of bias-corrected WALS estimates
			preserve 
			cap use `"`e(bcsimdata)'"', clear
			if _rc!=0 {
				noi di as err "{p 4 4 2}MC replications of bias-corrected wals estimates, "	/*
					*/ "as previously saved by {bf:`e(cmd)'}, do not exist. You need "	/*
					*/ "to refit the model.{p_end}" 
				exit _rc	
			}
			des , short varlist
			local wals_bc_list `r(varlist)'
			tempname wals_bc_repdata
			mata: `wals_bc_repdata'=st_data(., "`wals_bc_list'")
			restore
			local R=e(reps)
				
			* Confidence intervals for the conditional mean
			if "`type'"=="cinterval" {
				preserve
					keep if `touse'==1
					sort `e(ivar)' `tvar'
					tempname xf y_bar_i x_bar_i tV_i_s bcwp_rep 
					mata: `xf' = (st_data(., "`M_vars'"))'
					mata: `y_bar_i' = (st_data(., "`YBAR_i'"))'
					mata: `x_bar_i' = (st_data(., "`XBAR_list'"))'
					if "`e(model)'"=="fe-iid" {
						mata: `tV_i_s'=(`y_bar_i':*J(`R',`N',1))- `wals_bc_repdata'* `x_bar_i'
					}
					else {
						if e(Tcon)==1 {
							tempname ratio_m
							mata: `ratio_m'=st_numscalar("`ratio'")
							mata: `tV_i_s'=`ratio_m':*((`y_bar_i':*J(`R',`N',1))- `wals_bc_repdata'* `x_bar_i')
						}
						else {
							tempname RATIO_m
							mata: `RATIO_m' = (st_data(., "`RATIO'"))'
							mata: `tV_i_s'=(`RATIO_m':*J(`R',`N',1)):*((`y_bar_i':*J(`R',`N',1))- `wals_bc_repdata'* `x_bar_i')
						}
					}
					mata: `bcwp_rep'=`wals_bc_repdata'*`xf'+ `tV_i_s'
					local lev_lb	=(100-`level')/2
					local lev_ub	=100-`lev_lb'
					clear
					set obs `R'
					tempname bcwp_bounds 
					mata: `bcwp_bounds'=J(`N',2,0)
					tempvar bcwp_Rep
					gen double `bcwp_Rep'=.
					forvalues j=1(1)`N'{
						mata: st_store(.,., `bcwp_rep'[.,`j'])
						_pctile `bcwp_Rep', percentiles(`lev_lb' `lev_ub')
						mata: `bcwp_bounds'[`j',1]=st_numscalar("r(r1)")
						mata: `bcwp_bounds'[`j',2]=st_numscalar("r(r2)")
					}
				restore
				sort `e(ivar)' `tvar'
				tempvar bcwp_l bcwp_u
				gen double `bcwp_l'=. 
				gen double `bcwp_u'=.
				mata: st_store(.,st_varindex(tokens("`bcwp_l' `bcwp_u'")), st_varindex("`touse'"), `bcwp_bounds')
				noi gen `vtyp' `varn_lb' = `bcwp_l' if `touse'
				label var `varn_lb' "Lower bound `level'% confidence interval for conditional mean"
				noi gen `vtyp' `varn_ub' = `bcwp_u' if `touse'
				label var `varn_ub' "Upper bound `level'% confidence interval for conditional mean"
				exit
			}
			* Prediction intervals
			else if "`type'"=="pinterval" {
				tempname sig sig_M
				scalar `sig'=e(sigma)
				preserve
					keep if `touse'==1
					sort `e(ivar)' `tvar'
					tempname xf y_bar_i x_bar_i tV_i_s bcwp_rep 
					mata: `xf' = (st_data(., "`M_vars'"))'
					mata: `y_bar_i' = (st_data(., "`YBAR_i'"))'
					mata: `x_bar_i' = (st_data(., "`XBAR_list'"))'
					if "`e(model)'"=="fe-iid" {
						mata: `tV_i_s'=(`y_bar_i':*J(`R',`N',1))- `wals_bc_repdata'* `x_bar_i'
					}
					else {
						if e(Tcon)==1 {
							tempname ratio_m
							mata: `ratio_m'=st_numscalar("`ratio'")
							mata: `tV_i_s'=`ratio_m':*((`y_bar_i':*J(`R',`N',1))- `wals_bc_repdata'* `x_bar_i')
						}
						else {
							tempname RATIO_m
							mata: `RATIO_m' = (st_data(., "`RATIO'"))'
							mata: `tV_i_s'=(`RATIO_m':*J(`R',`N',1)):*((`y_bar_i':*J(`R',`N',1))- `wals_bc_repdata'* `x_bar_i')
						}
					}
					mata: `sig_M'=st_numscalar("`sig'")
					if "`rseed'"=="" {
						mata: rseed(strtoreal(st_local("rseed")))
					}
					mata: `bcwp_rep'=`wals_bc_repdata'*`xf'+ `tV_i_s' +`sig_M'*rnormal(`R',`N',0,1) 
					local lev_lb	=(100-`level')/2
					local lev_ub	=100-`lev_lb'
					clear
					set obs `R'
					tempname bcwp_bounds 
					mata: `bcwp_bounds'=J(`N',2,0)
					tempvar bcwp_Rep
					gen double `bcwp_Rep'=.
					forvalues j=1(1)`N'{
						mata: st_store(.,., `bcwp_rep'[.,`j'])
						_pctile `bcwp_Rep', percentiles(`lev_lb' `lev_ub')
						mata: `bcwp_bounds'[`j',1]=st_numscalar("r(r1)")
						mata: `bcwp_bounds'[`j',2]=st_numscalar("r(r2)")
					}
				restore
				sort `e(ivar)' `tvar'
				tempvar bcwp_l bcwp_u
				gen double `bcwp_l'=. 
				gen double `bcwp_u'=.
				mata: st_store(.,st_varindex(tokens("`bcwp_l' `bcwp_u'")), st_varindex("`touse'"), `bcwp_bounds')
				noi gen `vtyp' `varn_lb' = `bcwp_l' if `touse'
				label var `varn_lb' "Lower bound `level'% prediction interval"
				noi gen `vtyp' `varn_ub' = `bcwp_u' if `touse'
				label var `varn_ub' "Upper bound `level'% prediction interval"
				exit
			}
		}
	}
	error 198
end
*-------------------------------------------------------------------------------------------------

*! version 1.1
*! predictions of ar1wals
*! De Luca Giuseppe, Jan R. Magnus (2024/10/05)
program define ar1wals_p, sort
	if "`e(cmd)'" != "ar1wals" {
		noi di as error  "ar1wals was not the last estimation command"
		exit 301
	}
	version 14.0
	syntax [anything] [if] [in] [,xb biasp stdp rmsep bcp stdbcp Residuals BCResiduals wp bcwp CINTerval PINTerval LEVel(cilevel) rseed(numlist max=1 >0 integer)]
	if "`cinterval'"=="" & "`pinterval'"=="" {	
		
		* Parse options of the first syntax
		local myopts "biasp stdp rmsep bcp stdbcp Residuals BCResiduals wp bcwp"
		noi _pred_se "`myopts'" `0'
		if `s(done)' { 
			exit
		}
		local vtyp `s(typ)'
		local varn `s(varn)'
		local 0    `"`s(rest)'"'
		syntax [if] [in] [, `myopts']
		local type "`biasp'`stdp'`rmsep'`bcp'`stdbcp'`residuals'`bcresiduals'`wp'`bcwp'"
		
		* default: xb assumed 
		if ("`type'"=="") {
			di in gr "(option xb assumed; fitted values)"
			tempvar lp
			qui _predict double `lp' `if' `in', xb
			gen `vtyp' `varn' = `lp' `if' `in'
			label var `varn' "Linear prediction"
			exit
		}
		else if "`type'"=="biasp" {
			qui {
				tempname bias 
				matrix `bias'=e(bias)
				local Xlist: coln e(b)
				local ctype "`e(constype)'"
				tempvar bp
				if "`ctype'"!="noconstant" {
					local k:colsof e(bias)
					gen double `bp'=`bias'[1,`k'] 	`if' `in'
				}
				else {
					gen double `bp'=0				`if' `in'
				}    
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
		else if "`type'"=="stdp" {
			tempvar sep
			qui _predict double `sep' `if' `in', stdp
			gen `vtyp' `varn' = `sep' `if' `in'
			label var `varn' "SE of the prediction"
			exit
		}
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
		else if "`type'"=="stdbcp" {
			qui {
				local Xlist "`e(focvars)' `e(auxvars)'"				
				if "`e(constype)'"!="noconstant" {
					local Xlist: subinstr local Xlist "_cons" ""	
					local Xlist "`Xlist' _cons"	
				}
				preserve 
				tempname Vs
				use `"`e(bcsimdata)'"', clear
				corr, cov
				matrix `Vs'=r(C)
				restore
				if "`e(constype)'"!="noconstant" {
					tempvar const
					gen double `const'=1		`if' `in'
					local cons _cons
					local Xlist2: list Xlist - cons
					local Xlist2 "`Xlist2' `const'"
				}
				else {
					local Xlist2 "`Xlist'" 
				}
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
		else if "`type'"=="residuals" {
			tempvar lp
			qui _predict double `lp' `if' `in', xb
			gen `vtyp' `varn' = `e(depvar)'-`lp' `if' `in'
			label var `varn' "Residuals"
			exit
		}
		else if "`type'"=="bcresiduals" {
			tempvar bcp 
			qui predict double `bcp' `if' `in', bcp
			gen `vtyp' `varn' = `e(depvar)'-`bcp' `if' `in'
			label var `varn' "Bias-corrected residuals"
			exit
		}
		else if "`type'"=="wp" {
			qui {
				* Time variable and estimation sample
				tsset 
				local tvar "`r(timevar)'"
				tempvar esample  
				gen byte `esample'=e(sample)
				
				* List of non-collinear regressors 
				local Xlist "`e(focvars)' `e(auxvars)'"				
				if "`e(constype)'"!="noconstant" {
					local Xlist: subinstr local Xlist "_cons" ""	
				}
		
				* Mark sample (out-sample only)
				marksample touse
				markout `touse' `Xlist'
				replace `touse'=0 if `esample'==1

				* Last observation in the estimation sample 
				tempvar T_last  
				gsort -`esample' `tvar'
				gen byte `T_last'=`esample'==1 & _n==e(N)

				* Linear predictor, residual and powers for future periods
				tempvar lp wres pow 
				_predict double `lp' 	if `touse'==1, xb
				predict double `wres' 	if `T_last'==1, residuals
				tempname wresT_last
				sum `wres' 
				scalar `wresT_last'=r(mean)
				sum `tvar' if `T_last'==1
				gen double `pow'=`tvar'-r(mean) 	if `touse'==1
			}
			gen `vtyp' `varn' = `lp'+ (e(rho)^`pow')*`wresT_last' if `touse'==1
			label var `varn' "WALS prediction"
			exit
		}
		else if "`type'"=="bcwp" {
			qui {
				* Time variable and estimation sample
				tsset 
				local tvar "`r(timevar)'"
				tempvar esample bcp bcwres pow 
				gen byte `esample'=e(sample)

				* List of non-collinear regressors 
				local Xlist "`e(focvars)' `e(auxvars)'"				
				if "`e(constype)'"!="noconstant" {
					local Xlist: subinstr local Xlist "_cons" ""	
				}
		
				* Mark sample (out-sample only)
				marksample touse
				markout `touse' `Xlist'
				replace `touse'=0 if `esample'==1

				* Last observation in the estimation sample 
				tempvar T_last  
				gsort -`esample' `tvar'
				gen byte `T_last'=`esample'==1 & _n==e(N)

				* bc linear predictor, bc residual and powers for future periods
				tempvar bcp bcwres pow 
				predict double `bcp' 	if `touse'==1, bcp
				predict double `bcwres' if `T_last'==1, bcresiduals
				sum `tvar' if `T_last'==1
				gen double `pow'=`tvar'-r(mean) 	if `touse'==1
				tempname bcwresT_last
				sum `bcwres' 
				scalar `bcwresT_last'=r(mean)
			}
			gen `vtyp' `varn' = `bcp'+ (e(rho)^`pow')*`bcwresT_last' if `touse'==1
			label var `varn' "Bias-corrected WALS prediction"
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
				noi error 198
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

			* Time variable and estimation sample
			tsset 
			local tvar "`r(timevar)'"
			tempvar esample T_last lp wres pow 
			gen byte `esample'=e(sample)

			* List of non-collinear regressors 
			local Xlist "`e(focvars)' `e(auxvars)'"				
			if "`e(constype)'"!="noconstant" {
				tempname const
				gen double `const'=1 							
				local Xlist: subinstr local Xlist "_cons" ""	
				local Xlist "`Xlist' `const'"	
			}

			* List of regressors  (as temp vars)
			fvrevar `Xlist'							
			local M_vars "`r(varlist)'"

			* Mark sample (out-sample only)
			marksample touse
			markout `touse' `Xlist'
			replace `touse'=0 if `esample'==1
			count if `touse'==1
			local N=r(N)
		
			* Last observation in the estimation sample 
			tempvar T_last  
			gsort -`esample' `tvar'
			gen byte `T_last'=`esample'==1 & _n==e(N)

			* Powers for future periods
			tempvar pow 
			sum `tvar' if `T_last'==1
			gen double `pow'=`tvar'-r(mean) 	if `touse'==1
		
			* MC replications of bc WALS estimates
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

			* vector of MC replications of bc residuals at time T 
			tempname yT xT resT_rep
			gsort -`esample' `tvar'
			mata: `yT' = st_data(`e(N)', "`e(depvar)'")
			mata: `xT' = st_data(`e(N)', "`M_vars'")
			mata: `resT_rep'=J(`R',1,`yT')-`wals_bc_repdata'*(`xT'')  

			* Confidence intervals for E(y_{T+s}|x_{T+s})
			if "`type'"=="cinterval" {
				preserve
					keep if `touse'==1
					tempname xf Mpow rho_pow bcwp_rep 
					mata: `xf' = (st_data(., "`M_vars'"))'
					mata: `Mpow' = (st_data(., "`pow'"))'
					mata: `rho_pow' = `e(rho)':^`Mpow' 
					mata: `bcwp_rep'=`wals_bc_repdata'*`xf'+`resT_rep'*`rho_pow'  
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
				tempvar bcwp_l bcwp_u
				gen double `bcwp_l'=. 
				gen double `bcwp_u'=.
				mata: st_store(.,st_varindex(tokens("`bcwp_l' `bcwp_u'")), st_varindex("`touse'"), `bcwp_bounds')
				noi gen `vtyp' `varn_lb' = `bcwp_l' 	if `touse'
				label var `varn_lb' "Lower bound `level'% confidence interval for the conditional mean"
				noi gen `vtyp' `varn_ub' = `bcwp_u' if `touse'
				label var `varn_ub' "Upper bound `level'% confidence interval for the conditional mean"				
				exit
			}
			* Prediction interval for Y_{T+s}|X_{T+s}
			else if "`type'"=="pinterval" {
				tempname sig
				scalar `sig'=e(sigma)
				preserve
					keep if `touse'==1
					tempname xf Mpow rho_pow bcwp_rep 
					mata: `xf' = (st_data(., "`M_vars'"))'
					mata: `Mpow' = (st_data(., "`pow'"))'
					mata: `rho_pow' = `e(rho)':^`Mpow' 
					mata: sig=st_numscalar("`sig'")
					if "`rseed'"=="" {
						mata: rseed(strtoreal(st_local("rseed")))
					}
					mata: `bcwp_rep'=`wals_bc_repdata'*`xf'+`resT_rep'*`rho_pow' +sig*rnormal(`R',`N',0,1)  
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
				tempvar bcwp_l bcwp_u
				gen double `bcwp_l'=. 
				gen double `bcwp_u'=.
				mata: st_store(.,st_varindex(tokens("`bcwp_l' `bcwp_u'")), st_varindex("`touse'"), `bcwp_bounds')
				noi gen `vtyp' `varn_lb' = `bcwp_l' 	if `touse'
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

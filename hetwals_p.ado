*! version 1.0
*! predictions for hetwals
*! De Luca Giuseppe, Jan R. Magnus (2024/10/05)
program define hetwals_p, sort
	if "`e(cmd)'" != "hetwals" {
		noi di as error  "hetwals was not the last estimation command"
		exit 301
	}
	version 14.0
	syntax [anything] [if] [in] [,xb biasp stdp rmsep bcp stdbcp Residuals BCResiduals sigma CINTerval PINTerval LEVel(cilevel) rseed(numlist max=1 >0 integer)]
	if "`cinterval'"=="" & "`pinterval'"=="" {	
		
		local myopts "biasp stdp rmsep bcp stdbcp Residuals BCResiduals sigma"
		noi _pred_se "`myopts'" `0'
		if `s(done)' { 
			exit
		}
		local vtyp `s(typ)'
		local varn `s(varn)'
		local 0    `"`s(rest)'"'
		syntax [if] [in] [, `myopts']
		local type "`biasp'`stdp'`rmsep'`bcp'`stdbcp'`residuals'`bcresiduals'`sigma'"
		
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
		else if "`type'"=="sigma" {
			qui {
				tempname hb 
				matrix `hb'=e(hb)
				local Vlist: coln e(hb)
				local q: colsof e(hb)
				local cc "_cons"
				local Vlist: list Vlist-cc 
				tempvar lnvar
				gen double `lnvar'=`hb'[1,`q'] 	`if' `in'
				local h 1
				foreach v of local Vlist {
					qui replace `lnvar'=`lnvar'+`hb'[1,`h']*`v'
					local h=`h'+1
				}
			}
			gen `vtyp' `varn' = exp(0.5 * `lnvar')
			label var `varn' "Estimated standard deviations"
			exit
		}
	}
	else {
		qui {
			* Parse options 
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

			* List of non-collinear covariates 
			local Xlist "`e(focvars)' `e(auxvars)'"				
			if "`e(constype)'"!="noconstant" {
				tempname const
				gen double `const'=1 							
				local Xlist: subinstr local Xlist "_cons" ""	
				local Xlist "`Xlist' `const'"	
			}
			
			* Mark sample
			marksample touse
			markout `touse' `Xlist'  
			count if `touse'==1
			local N=r(N)
			
			* List of focus and auxliary regressors in Mata (as temp vars)
			fvrevar `Xlist'								if `touse' 
			local M_vars "`r(varlist)'"
			
			* Dataset of MC replications 
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

			* Confidence intervals for conditional mean
			if "`type'"=="cinterval" {
				preserve
					keep if `touse'==1
					tempname xf xf_rep
					mata: `xf' = st_data(., "`M_vars'")
					mata: `xf_rep'=`wals_bc_repdata'*(`xf'')  
					local lev_lb	=(100-`level')/2
					local lev_ub	=100-`lev_lb'
					clear
					set obs `R'
					tempname XF_rep xf_bounds 
					mata: `xf_bounds'=J(`N',2,0)
					gen double `XF_rep'=.
					forvalues j=1(1)`N'{
						mata: st_store(.,., `xf_rep'[.,`j'])
						_pctile `XF_rep', percentiles(`lev_lb' `lev_ub')
						mata: `xf_bounds'[`j',1]=st_numscalar("r(r1)")
						mata: `xf_bounds'[`j',2]=st_numscalar("r(r2)")
					}
				restore
				tempname xf_l xf_u
				gen double `xf_l'=. 
				gen double `xf_u'=.
				mata: st_store(.,st_varindex(tokens("`xf_l' `xf_u'")), st_varindex("`touse'"), `xf_bounds')
				noi gen `vtyp' `varn_lb' = `xf_l' if `touse'
				label var `varn_lb' "Lower bound `level'% confidence interval for conditional mean"
				noi gen `vtyp' `varn_ub' = `xf_u' if `touse'
				label var `varn_ub' "Upper bound `level'% confidence interval for conditional mean"
				exit
			}
			
			* Prediction intervals
			else if "`type'"=="pinterval" {
				tempname hb 
				matrix `hb'=e(hb)
				local Vlist: coln e(hb)
				local q: colsof e(hb)
				local cc "_cons"
				local Vlist: list Vlist-cc 
				tempvar lnvar hat_sig
				gen double `lnvar'=`hb'[1,`q'] 			if `touse'
				local h 1
				foreach v of local Vlist {
					replace `lnvar'=`lnvar'+`hb'[1,`h']*`v'
					local h=`h'+1
				}
				gen double `hat_sig' = exp(0.5 * `lnvar')	if `touse'
				
				preserve
					keep if `touse'==1
					tempname sig_i xf xf_rep 
					mata: `sig_i'= st_data(., "`hat_sig'")	
					mata: `xf' = st_data(., "`M_vars'")
					if "`rseed'"=="" {
						mata: rseed(strtoreal(st_local("rseed")))
					}
					mata: `xf_rep'=`wals_bc_repdata'*(`xf'')+`sig_i'':*rnormal(`R',`N'',0,1)  
					local lev_lb	=(100-`level')/2
					local lev_ub	=100-`lev_lb'
					clear
					set obs `R'
					tempname XF_rep xf_bounds 
					mata: `xf_bounds'=J(`N',2,0)
					gen double `XF_rep'=.
					forvalues j=1(1)`N'{
						mata: st_store(.,., `xf_rep'[.,`j'])
						_pctile `XF_rep', percentiles(`lev_lb' `lev_ub')
						mata: `xf_bounds'[`j',1]=st_numscalar("r(r1)")
						mata: `xf_bounds'[`j',2]=st_numscalar("r(r2)")
					}
				restore
				tempname xf_l xf_u
				gen double `xf_l'=. 
				gen double `xf_u'=.
				mata: st_store(.,st_varindex(tokens("`xf_l' `xf_u'")), st_varindex("`touse'"), `xf_bounds')
				noi gen `vtyp' `varn_lb' = `xf_l' if `touse'
				label var `varn_lb' "Lower bound `level'% prediction interval"
				noi gen `vtyp' `varn_ub' = `xf_u' if `touse'
				label var `varn_ub' "Upper bound `level'% prediction interval"
				exit
			}
		}
	}
	error 198
end
*-------------------------------------------------------------------------------------------------



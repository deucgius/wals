*! version 1.1
*! predictions of wals
*! De Luca Giuseppe (2024/10/05)
program define wals_p, sort
	if "`e(cmd)'" != "wals" {
		noi di as error  "wals was not the last estimation command"
		exit 301
	}
	version 14.0
	syntax [anything] [if] [in] [,xb biasp stdp rmsep bcp stdbcp Residuals BCResiduals CINTerval PINTerval LEVel(cilevel) rseed(numlist max=1 >0 integer)]
	if "`cinterval'"=="" & "`pinterval'"=="" {	
		
		* Parse options of the first syntax
		local myopts "Residuals stdp biasp rmsep bcp BCResiduals stdbcp"
		noi _pred_se "`myopts'" `0'
		if `s(done)' { 
			exit
		}
		local vtyp `s(typ)'
		local varn `s(varn)'
		local 0    `"`s(rest)'"'
		syntax [if] [in] [, `myopts']
		local type "`biasp'`stdp'`rmsep'`residuals'`bcp'`bcresiduals'`stdbcp'"
		
		* default: xb assumed 
		if ("`type'"=="") {
			di in gr "(option xb assumed; fitted values)"
			tempvar lp
			qui _predict double `lp' `if' `in', xb
			gen `vtyp' `varn' = `lp' `if' `in'
			label var `varn' "Linear prediction"
			exit
		}
		else if "`type'"=="residuals" {
			tempvar lp
			qui _predict double `lp' `if' `in', xb
			gen `vtyp' `varn' = `e(depvar)'-`lp' `if' `in'
			label var `varn' "Residuals"
			exit
		}
		else if "`type'"=="biasp" {
			tempname bias 
			matrix `bias'=e(bias)
			local Xlist: coln e(b)
			local ctype "`e(constype)'"
			tempvar bp
			if "`ctype'"!="noconstant" {
				local k:colsof e(bias)
				qui gen double `bp'=`bias'[1,`k'] 	`if' `in'
			}
			else {
				qui gen double `bp'=0				`if' `in'
			}    
			local h 1
			foreach v of local Xlist {
				if "`v'"!="_cons" qui replace `bp'=`bp'+`bias'[1,`h']*`v'
				local h=`h'+1
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
			tempname bias 
			matrix `bias'=e(bias)
			local Xlist: coln e(b)
			local ctype "`e(constype)'"
			tempvar bp
			if "`ctype'"!="noconstant" {
				local k:colsof e(bias)
				qui gen double `bp'=`bias'[1,`k'] 	`if' `in'
			}
			else {
				qui gen double `bp'=0				`if' `in'
			}    
			local h 1
			foreach v of local Xlist {
				if "`v'"!="_cons" qui replace `bp'=`bp'+`bias'[1,`h']*`v'
				local h=`h'+1
			}
			tempvar sep
			qui _predict double `sep' `if' `in', stdp
			gen `vtyp' `varn' = sqrt(`bp'^2+`sep'^2) `if' `in'
			label var `varn' "RMSE of the prediction"
			exit
		}
		else if "`type'"=="bcp" {
			tempvar lp
			qui _predict double `lp' `if' `in', xb
			tempname bias 
			matrix `bias'=e(bias)
			local Xlist: coln e(b)
			local ctype "`e(constype)'"
			tempvar bp
			if "`ctype'"!="noconstant" {
				local k:colsof e(bias)
				qui gen double `bp'=`bias'[1,`k'] 	`if' `in'
			}
			else {
				qui gen double `bp'=0				`if' `in'
			}    
			local h 1
			foreach v of local Xlist {
				if "`v'"!="_cons" qui replace `bp'=`bp'+`bias'[1,`h']*`v'
				local h=`h'+1
			}
			gen `vtyp' `varn' = `lp'-`bp' `if' `in'
			label var `varn' "Bias-corrected prediction"
			exit
		}
		else if "`type'"=="bcresiduals" {
			tempvar lp
			qui _predict double `lp' `if' `in', xb
			tempname bias 
			matrix `bias'=e(bias)
			local Xlist: coln e(b)
			local ctype "`e(constype)'"
			tempvar bp
			if "`ctype'"!="noconstant" {
				local k:colsof e(bias)
				qui gen double `bp'=`bias'[1,`k'] 	`if' `in'
			}
			else {
				qui gen double `bp'=0				`if' `in'
			}    
			local h 1
			foreach v of local Xlist {
				if "`v'"!="_cons" qui replace `bp'=`bp'+`bias'[1,`h']*`v'
				local h=`h'+1
			}
			gen `vtyp' `varn' = `e(depvar)'-`lp'+`bp' `if' `in'
			label var `varn' "Bias-corrected residuals"
			exit
		}
		else if "`type'"=="stdbcp" {
			local Xlist "`e(focvars)' `e(auxvars)'"				
			if "`e(constype)'"!="noconstant" {
				local Xlist: subinstr local Xlist "_cons" ""	
				local Xlist "`Xlist' _cons"	
			}
			preserve 
			tempname Vs
			use `"`e(bcsimdata)'"', clear
			qui corr, cov
			matrix `Vs'=r(C)
			restore
			if "`e(constype)'"!="noconstant" {
				tempvar const
				qui gen double `const'=1		`if' `in'
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
			qui count if `touse'==1
			local N=r(N)
			tempvar bcs
			qui gen double `bcs'=.	
			tempname xj xjVsxjt
			forvalues j=1(1)`N' {
				mkmat `Xlist2' if _n==`j', matrix(`xj')
				matrix `xjVsxjt'=`xj' * `Vs' * `xj''
				qui replace `bcs'=sqrt(`xjVsxjt'[1,1]) if _n==`j'
			}
			gen `vtyp' `varn' = `bcs' if `touse'
			label var `varn' "S.E. of bias-corrected prediction"
			exit
		}
	}
	else {
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
			di as err "only one of {bf:cinterval} or {bf:pinterval} is allowed"
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
			qui gen double `const'=1 							
			local Xlist: subinstr local Xlist "_cons" ""	
			local Xlist "`Xlist' `const'"	
		}
		
		* Mark sample
		marksample touse
		markout `touse' `Xlist'  
		qui count if `touse'==1
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
		qui des , short varlist
		local wals_bc_list `r(varlist)'
		tempname wals_bc_repdata
		mata: `wals_bc_repdata'=st_data(., "`wals_bc_list'")
		restore
		local R=e(reps)

		* Confidence intervals for E(y|x)
		if "`type'"=="cinterval" {
			tempname sig
			scalar `sig'=e(sigma)
			preserve
				qui keep if `touse'==1
				tempname xf xf_rep
				qui mata: `xf' = st_data(., "`M_vars'")
				qui mata: `xf_rep'=`wals_bc_repdata'*(`xf'')  
				local lev_lb	=(100-`level')/2
				local lev_ub	=100-`lev_lb'
				clear
				qui set obs `R'
				tempname XF_rep xf_bounds 
				mata: `xf_bounds'=J(`N',2,0)
				qui gen double `XF_rep'=.
				forvalues j=1(1)`N'{
					qui mata: st_store(.,., `xf_rep'[.,`j'])
					_pctile `XF_rep', percentiles(`lev_lb' `lev_ub')
					qui mata: `xf_bounds'[`j',1]=st_numscalar("r(r1)")
					qui mata: `xf_bounds'[`j',2]=st_numscalar("r(r2)")
				}
			restore
			tempname xf_l xf_u
			qui gen double `xf_l'=. 
			qui gen double `xf_u'=.
			qui mata: st_store(.,st_varindex(tokens("`xf_l' `xf_u'")), st_varindex("`touse'"), `xf_bounds')
			gen `vtyp' `varn_lb' = `xf_l' if `touse'
			label var `varn_lb' "Lower bound `level'% prediction interval"
			qui gen `vtyp' `varn_ub' = `xf_u' if `touse'
			label var `varn_ub' "Upper bound `level'% prediction interval"
			exit
		}
		
		* Prediction intervals for y 
		else if "`type'"=="pinterval" {
			tempname sig
			scalar `sig'=e(sigma)
			preserve
				qui keep if `touse'==1
				tempname xf xf_rep
				qui mata: sig=st_numscalar("`sig'")
				qui mata: `xf' = st_data(., "`M_vars'")
				if "`rseed'"=="" {
					qui mata: rseed(strtoreal(st_local("rseed")))
				}
				qui mata: `xf_rep'=`wals_bc_repdata'*(`xf'')+sig*rnormal(`R',`N'',0,1)  
				local lev_lb	=(100-`level')/2
				local lev_ub	=100-`lev_lb'
				clear
				qui set obs `R'
				tempname XF_rep xf_bounds 
				mata: `xf_bounds'=J(`N',2,0)
				qui gen double `XF_rep'=.
				forvalues j=1(1)`N'{
					qui mata: st_store(.,., `xf_rep'[.,`j'])
					_pctile `XF_rep', percentiles(`lev_lb' `lev_ub')
					qui mata: `xf_bounds'[`j',1]=st_numscalar("r(r1)")
					qui mata: `xf_bounds'[`j',2]=st_numscalar("r(r2)")
				}
			restore
			tempname xf_l xf_u
			qui gen double `xf_l'=. 
			qui gen double `xf_u'=.
			qui mata: st_store(.,st_varindex(tokens("`xf_l' `xf_u'")), st_varindex("`touse'"), `xf_bounds')
			gen `vtyp' `varn_lb' = `xf_l' if `touse'
			label var `varn_lb' "Lower bound `level'% prediction interval"
			qui gen `vtyp' `varn_ub' = `xf_u' if `touse'
			label var `varn_ub' "Upper bound `level'% prediction interval"
			exit
		}
	}
	error 198
end
*-------------------------------------------------------------------------------------------------

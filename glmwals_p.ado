*! version 1.0
*! predictions for glmwals
*! De Luca Giuseppe, Jan R. Magnus (2024/10/05)
program define glmwals_p, sort
	if "`e(cmd)'" != "glmwals" {
		noi di as error  "glmwals was not the last estimation command"
		exit 301
	}
	version 14.0
	syntax [anything] [if] [in] [, mu xb biasp stdp rmsep bcp stdbcp Residuals]
	local myopts "mu biasp stdp rmsep bcp stdbcp Residuals"
	noi _pred_se "`myopts'" `0'
	if `s(done)' { 
		exit
	}
	local vtyp `s(typ)'
	local varn `s(varn)'
	local 0    `"`s(rest)'"'
	syntax [if] [in] [, `myopts']
	local type "`mu'`biasp'`stdp'`rmsep'`bcp'`stdbcp'`residuals'"
	
	* default: mu assumed 
	if "`type'"==""|"`type'"=="mu" {
		if "`type'"=="" {
			di in gr "(option mu assumed; predicted mean `e(depvar)')"
		}
		qui{ 
			tempvar lp mu
			_predict double `lp' `if' `in', xb `offset'
			global GLMWALS_m	`"`e(btrials)'"'		/* Binomial denominator */
			global GLMWALS_p	`"`e(power)'"'			/* Power   			 	*/
			`e(link)' 1 `lp' `mu' 		
			macro drop GLMWALS_*
		}
		gen `vtyp' `varn' = `mu' `if' `in'
		label var `varn' "Predicted mean `e(depvar)'"
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
		tempvar mu
		qui predict double `mu' `if' `in', mu
		gen `vtyp' `varn' = `e(depvar)'-`mu' `if' `in'
		label var `varn' "Residuals"
		exit
	}
	error 198
end
*-------------------------------------------------------------------------------------------------

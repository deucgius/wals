//-----------------------------------------------------------------------------------------------------------------------	
// This program is a modified version of glim_v2.ado (version 1.2.0, 19feb2019), 
//		which defines the Binomial variance/family 
//-----------------------------------------------------------------------------------------------------------------------	
program define glmwals_v2
	version 14.0
	args todo eta mu return 
	if `todo' == -1 {													/* Title 		*/
		local y     "${GLMWALS_y}"
		local m     "${GLMWALS_m}"
		local touse "`eta'"
		capture assert `m'>0 if `touse'         
		if _rc {
			di as err `"`m' has nonpositive values"'
			exit 499
		}
		capture assert `m'==int(`m') if `touse'
		if _rc {
			capture confirm number `m'
			if _rc {
				di as txt `"note: `m' has noninteger values"'
			}
			else {
				di as txt `"note: `m' is a noninteger"'
			}
		}
		capture assert `y'>=0 if `touse'
		if _rc {
			di as err `"dependent variable `y' has negative values"'
			exit 499
		}
		capture assert `y'<=`m' if `touse' & `y' != .
		if _rc {
			di as err `"`y' > `m' in some cases"'
			exit 499
		}
		capture assert `y'==int(`y') if `touse'
		if _rc {
		 	di as txt `"note: `y' has noninteger values"'	
		}
		cap assert `y'== 0 | `y' == 1
		if "`m'"!="1" | _rc {			
			global GLMWALS_vt "Binomial"
			global GLMWALS_vf "u*(1-u/`m')"
		}
		else {
			global GLMWALS_vt "Bernoulli"
			global GLMWALS_vf "u*(1-u)"
		}
		global GLMWALS_mu "glmwals_mu 0 `m'"
		exit
	}
	if `todo' == 1 {													/* V(mu) 			*/
		gen double `return' =  `mu'*(1-`mu'/${GLMWALS_m})
		exit 
	}
	noi di as err "Unknown call to glmwals variance function"
	error 198
end
//-----------------------------------------------------------------------------------------------------------------------	

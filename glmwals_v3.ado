//-----------------------------------------------------------------------------------------------------------------------	
// This program is a modified version of glim_v3.ado (version 1.1.0, 02apr2009), 
//		which defines the Poisson variance/family 
//-----------------------------------------------------------------------------------------------------------------------	
program define glmwals_v3
	version 14.0
	args todo eta mu return
	if `todo' == -1 {													/* Title 			*/
		local y     "${GLMWALS_y}"
		local touse "`eta'"
		capture assert `y'>=0 if `touse'
		if _rc {
			di as err `"dependent variable `y' has negative values"'
			exit 499
		}
		capture assert `y'==int(`y') if `touse'
		if _rc {
			di as txt `"note: `y' has noninteger values"'
		}
		global GLMWALS_vt "Poisson"
		global GLMWALS_vf "u"
		global GLMWALS_mu "glmwals_mu 0 ."	
		exit
	}
	if `todo' == 1 {													/* V(mu) 			*/
		gen double `return' =  `mu'
		exit 
	}
	noi di as err "Unknown call to glmwals variance function"
	error 198
end
//-----------------------------------------------------------------------------------------------------------------------	


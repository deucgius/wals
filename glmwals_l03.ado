//-----------------------------------------------------------------------------------------------------------------------	
// This program is a modified version of glim_l03.ado (version 1.0.1, 28apr2000), 
//		which defines the "log" link
//-----------------------------------------------------------------------------------------------------------------------	
program define glmwals_l03				
	version 14.0
	args todo eta mu return
	if `todo' == -1 {                       							/* Title */
		global GLMWALS_lt "Log"
		if "${GLMWALS_m}" == "1" {
			global GLMWALS_lf "ln(u)"
		}
		else {
			global GLMWALS_lf "ln(u/${GLMWALS_m})"
		}
		exit
	}
	if `todo' == 0 {													/* eta = g(mu) */
		gen double `eta' = ln(`mu'/${GLMWALS_m})
		exit 
	}
	if `todo' == 1 {													/* mu = g^-1(eta) */
		gen double `mu' = ${GLMWALS_m}*exp(`eta')
		exit 
	}
	if `todo' == 2 {													/* (d mu)/(d eta) */
		gen double `return' = `mu'
		exit 
	}
	noi di as err "Unknown call to glmwals link function"
	exit 198
end
//-----------------------------------------------------------------------------------------------------------------------


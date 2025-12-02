//-----------------------------------------------------------------------------------------------------------------------	
// This program is a modified version of glim_l05.ado (version 1.1.0, 19feb2019), 
//		which defines the "log complement" link 
//-----------------------------------------------------------------------------------------------------------------------	
program define glmwals_l05				
	version 14.0
	args todo eta mu return
	if `todo' == -1 {                       							/* Title */
		global SGLM_lt "Log complement"
		if "${GLMWALS_m}" == "1" {
			global GLMWALS_lf "ln(1-u)"
		}
		else {
			global GLMWALS_lf "ln(1-u/${GLMWALS_m})"
		}
		exit
	}
	if `todo' == 0 {													/* eta = g(mu) */
		gen double `eta' = ln1m(`mu'/${GLMWALS_m})
		exit 
	}
	if `todo' == 1 {													/* mu = g^-1(eta) */
		gen double `mu' = ${GLMWALS_m}*(-expm1(`eta'))
		exit 
	}
	if `todo' == 2 {													/* (d mu)/(d eta) */
		gen double `return' = `mu'-${GLMWALS_m}
		exit 
	}
	noi di as err "Unknown call to glmwals link function"
	exit 198
end
//-----------------------------------------------------------------------------------------------------------------------


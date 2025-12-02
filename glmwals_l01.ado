//-----------------------------------------------------------------------------------------------------------------------	
// This program is a modified version of glim_l01.ado (version 1.0.1, 28apr2000), 
//		which defines the "identity" link
//-----------------------------------------------------------------------------------------------------------------------	
program define glmwals_l01							
	version 14.0
	args todo eta mu return
	if `todo' == -1 {                      								/* Title */
		global GLMWALS_lt "Identity"
		if "${GLMWALS_m}" == "1" {	
			global GLMWALS_lf "u"
		}
		else {
			global GLMWALS_lf "u/${GLMWALS_m}"
		}
		exit
	}
	if `todo' == 0 {													/* eta = g(mu) */
		gen double `eta' = `mu'/${GLMWALS_m}
		exit 
	}
	if `todo' == 1 {													/* mu = g^-1(eta) */
		gen double `mu' = `eta'*${GLMWALS_m}
		exit 
	}
	if `todo' == 2 {													/* (d mu)/(d eta) */
		gen double `return' = ${GLMWALS_m}
		exit 
	}
	noi di as err "Unknown call to glmwals link function"
	exit 198
end
//-----------------------------------------------------------------------------------------------------------------------


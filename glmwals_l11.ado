//-----------------------------------------------------------------------------------------------------------------------	
// This program is a modified version of glim_l11.ado (version 1.0.1, 28apr2000), 
//		which defines the "Power(a)" link 
//-----------------------------------------------------------------------------------------------------------------------	
program define glmwals_l11			
	version 14.0
	args todo eta mu return
	if `todo' == -1 {                       							/* Title */
		global GLMWALS_lt "Power"
		if "${GLMWALS_m}" == "1" {
			global GLMWALS_lf "u^(${GLMWALS_p})"
		}
		else {
			global GLMWALS_lf "(u/${GLMWALS_m})^(${GLMWALS_p})"
		}
		exit
	}
	if `todo' == 0 {													/* eta = g(mu) */
		gen double `eta' = (`mu'/${GLMWALS_m})^(${GLMWALS_p})
		exit 
	}
	if `todo' == 1 {													/* mu = g^-1(eta) */
		gen double `mu' = ${GLMWALS_m}*`eta'^(1/${GLMWALS_p})
		exit 
	}
	if `todo' == 2 {													/* (d mu)/(d eta) */
		gen double `return' = `mu'^(1-${GLMWALS_p})*${GLMWALS_m}^${GLMWALS_p}/${GLMWALS_p}
		exit 
	}
	noi di as err "Unknown call to glmwals link function"
	exit 198
end
//-----------------------------------------------------------------------------------------------------------------------	


//-----------------------------------------------------------------------------------------------------------------------	
// This program is a modified version of glim_l12.ado (version 1.0.1, 28apr2000), 
//		which defines the "Odds power(a)" link 
//-----------------------------------------------------------------------------------------------------------------------	
program define glmwals_l12				
	version 14.0
	args todo eta mu return
	if `todo' == -1 {                       							/* Title */
		global GLMWALS_lt "Odds power"
		if "${GLMWALS_m}" == "1" {
			global GLMWALS_lf "((u/(1-u))^(${GLMWALS_p})-1)/(${GLMWALS_p})"
		}
		else {
			global GLMWALS_lf "((u/(${GLMWALS_m}-u))^(${GLMWALS_p})-1)/(${GLMWALS_p})"
		}
        exit
	}
	if `todo' == 0 {													/* eta = g(mu) */
		gen double `eta' = ((`mu'/(${GLMWALS_m}-`mu'))^(${GLMWALS_p})-1)/${GLMWALS_p}
		exit 
	}
	if `todo' == 1 {													/* mu = g^-1(eta) */
		gen double `mu' = ${GLMWALS_m}*(1+${GLMWALS_p}*`eta')^(1/${GLMWALS_p}) / /*
				*/ (1 + (1+${GLMWALS_p}*`eta')^(1/${GLMWALS_p})) 
		exit 
	}
	if `todo' == 2 {													/* (d mu)/(d eta) */
		gen double `return' = ${GLMWALS_m}*(`mu'/${GLMWALS_m})^(1-${GLMWALS_p}) * /*
			*/ (1-`mu'/${GLMWALS_m})^(1+${GLMWALS_p})
		exit 
	}
	noi di as err "Unknown call to glmwals link function"
	exit 198
end
//-----------------------------------------------------------------------------------------------------------------------	

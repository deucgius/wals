//-----------------------------------------------------------------------------------------------------------------------	
// This program is a modified version of glim_l08.ado (version 1.0.1, 28apr2000), 
//		which defines the "Probit" link 
//-----------------------------------------------------------------------------------------------------------------------	
program define glmwals_l08				
	version 14.0
	args todo eta mu return
	if `todo' == -1 {                       							/* Title */
		global GLMWALS_lt "Probit"
		if "${GLMWALS_m}" == "1" {	
			global GLMWALS_lf "invnormal(u)"
		}
		else {
			global GLMWALS_lf "invnormal(u/${GLMWALS_m})"
		}
		exit
	}
	if `todo' == 0 {													/* eta = g(mu) */
		gen double `eta' = invnormal(`mu'/${GLMWALS_m})
		exit 
	}
	if `todo' == 1 {													/* mu = g^-1(eta) */
		gen double `mu' = ${GLMWALS_m}*normal(`eta')
		exit 
	}
	if `todo' == 2 {													/* (d mu)/(d eta) */
		gen double `return' = ${GLMWALS_m}*normalden(invnormal(`mu'/${GLMWALS_m}))
		exit 
	}
	noi di as err "Unknown call to glmwals link function"
	exit 198
end
//-----------------------------------------------------------------------------------------------------------------------	

//-----------------------------------------------------------------------------------------------------------------------	
// This program is a modified version of glim_l10.ado (version 1.0.1, 28apr2000), 
//		which defines the "Power(-2)" link 
//-----------------------------------------------------------------------------------------------------------------------	
program define glmwals_l10			
	version 14.0
	args todo eta mu return
	if `todo' == -1 {                       							/* Title */
		global GLMWALS_lt "Power(-2)"
		global GLMWALS_lf "1/(u^2)"
		exit
	}
	if `todo' == 0 {													/* eta = g(mu) */
		gen double `eta' = (${GLMWALS_m}/`mu')*(${GLMWALS_m}/`mu')
		exit 
	}
	if `todo' == 1 {													/* mu = g^-1(eta) */
		gen double `mu' = ${GLMWALS_m}/sqrt(`eta')
		exit 
	}
	if `todo' == 2 {													/* (d mu)/(d eta) */
		gen double `return' = -(`mu'*`mu'*`mu')/(2*${GLMWALS_m}*${GLMWALS_m})
		exit 
	}
	noi di as err "Unknown call to glmwals link function"
	exit 198
end
//-----------------------------------------------------------------------------------------------------------------------	

//-----------------------------------------------------------------------------------------------------------------------	
// This program is a modified version of glim_l09.ado (version 1.0.1, 28apr2000), 
//		which defines the "Reciprocal" link 
//-----------------------------------------------------------------------------------------------------------------------	
program define glmwals_l09				
	version 14.0
	args todo eta mu return
	if `todo' == -1 {                      								/* Title */
		global GLMWALS_lt "Reciprocal"
		if "${GLMWALS_m}" == "1" {
			global GLMWALS_lf "1/u"
		}
		else {
			global GLMWALS_lf "${GLMWALS_m}/u"
		}
        exit
	}
	if `todo' == 0 {													/* eta = g(mu) */
		gen double `eta' = ${GLMWALS_m}/`mu'
		exit 
	}
	if `todo' == 1 {													/* mu = g^-1(eta) */
		gen double `mu' = ${GLMWALS_m}/`eta'
		exit 
	}
	if `todo' == 2 {													/* (d mu)/(d eta) */
		gen double `return' = -`mu'*`mu'/${GLMWALS_m}
		exit 
	}
	noi di as err "Unknown call to glmwals link function"
	exit 198
end
//-----------------------------------------------------------------------------------------------------------------------	


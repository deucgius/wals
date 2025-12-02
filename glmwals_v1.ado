//-----------------------------------------------------------------------------------------------------------------------	
// This program is a modified version of glim_v1.ado (version 1.1.0, 02apr2009), 
//		which defines the Gaussian variance/family 
//-----------------------------------------------------------------------------------------------------------------------	
program define glmwals_v1
	version 14.0
	args todo eta mu return
	if `todo' == -1 {													/* Title 			*/
		global GLMWALS_vt "Gaussian"
		global GLMWALS_vf "1"
		global GLMWALS_mu "*"		
		exit
	}
	if `todo' == 1 {													/* V(mu) 			*/
		scalar `return' = 1 
		exit 
	}
	noi di as err "Unknown call to glmwals variance function"
	error 198
end
//-----------------------------------------------------------------------------------------------------------------------	

//-----------------------------------------------------------------------------------------------------------------------	
// This program is a modified version of glim_mu.ado (version 1.0.0, 28apr2000), 
//		which shifts (if needed) the mean outcome into its range 
// 		when this is not ensured by the link 
//-----------------------------------------------------------------------------------------------------------------------	
program define glmwals_mu
	version 14.0
	args min max mu
	tempname eps
	scalar `eps' = 1e-8
	replace `mu' = `min'+`eps' if `mu'<=`min'+`eps'
	replace `mu' = `max'-`eps' if `mu'>=`max'-`eps'
end
//-----------------------------------------------------------------------------------------------------------------------	

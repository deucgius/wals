//-----------------------------------------------------------------------------------------------------------------------	
// This program compute confidence intervals using replicated of bias-corrected wals estimates
//-----------------------------------------------------------------------------------------------------------------------	
program define walsci							
qui{
	version 14.0
	args data lev return

	preserve 
	capture use `"`data'"', clear
	if _rc!=0 {
		noi di as err "Simulation data not fount"
		exit 198
	}
	capture assert `lev'>=10 & `lev'<=99.99
	if _rc!=0 {
		noi di as res "`lev'"
		noi di as err "level() must be between 10 and 99.99 inclusive"
		exit 198
	}
	local K=`c(k)'
	ds
	unab varlist : _all
	tokenize "`varlist'"
	local bc_var: subinstr local 1 "_1" ""
	local vars ""
	local lev_lb	=(100-`lev')/2
	local lev_ub	=100-`lev_lb'
	matrix `return'=J(2,`K',.)
	forvalues h=1(1)`K'{
		_pctile `bc_var'_`h' , percentiles(`lev_lb' `lev_ub')
		matrix `return'[1,`h']=r(r1)
		matrix `return'[2,`h']=r(r2)
		local vv: variable label `bc_var'_`h'
		local vars "`vars' `vv'"
	}
	matrix rown `return' = ll ul
	matrix coln `return' = `vars' 
	restore
}
end
//-----------------------------------------------------------------------------------------------------------------------

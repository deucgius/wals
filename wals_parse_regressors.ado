//-----------------------------------------------------------------------------------------------------------------------	
// wals_parse_regressors: parse focus and auxiliary regressors
//-----------------------------------------------------------------------------------------------------------------------	
program wals_parse_regressors, sclass
	syntax [anything]
	sreturn clear
	local foclist0 ""
	local auxlist0 ""
	local fvindvars `anything'
	while `"`fvindvars'"' != "" {
		gettoken varlist1 fvindvars: fvindvars, bind
		gettoken varlist : varlist1, bind match(paren)
		if "`paren'" == "(" {
			local 0 `varlist'	
			cap syntax [anything]
			if _rc!=0 {
				di as err "focus regressors must be specified in parentheses without options"
				exit 198
			}
			if "`anything'" == "" {
				di as err "varlist may not be empty: () not allowed"
				exit 198
			}
			syntax varlist(numeric fv ts)
			local foclist0 `foclist0' `varlist'
		}
		else  {
			local 0 `varlist'	
			syntax varlist(numeric fv ts)
			local auxlist0 `auxlist0' `varlist'
		}
	}
	sreturn local foc  = `"`foclist0'"'
	sreturn local aux  = `"`auxlist0'"'
end
//-----------------------------------------------------------------------------------------------------------------------	

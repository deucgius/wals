* Estimating linear combination of parameters after WALS regression
*! Version 1.0
* Date			: 2024/10/05
* Authors		: De Luca Giuseppe, Jan R. Magnus
*-------------------------------------------------------------------------------------------------
* lcwals
*-------------------------------------------------------------------------------------------------
program define lcwals, rclass 
version 14.0, missing
qui {
	if 		"`e(cmd)'"!="wals"		&	/*
		*/	"`e(cmd)'"!="hetwals"	&	/*
		*/	"`e(cmd)'"!="ar1wals"	&	/*
		*/	"`e(cmd)'"!="xtwals"	&	/*
		*/	"`e(cmd)'"!="glmwals" {
		noi di as error  "wals was not the last estimation command"
		exit 301
	}

	* Parse formula for the linear combination
	gettoken token 0 : 0, parse(",= ") bind
	while `"`token'"'!="" & `"`token'"'!="," {
		if `"`token'"' == "=" {
			noi di in red _quote "=" _quote /*
			*/ " not allowed in expression"
			exit 198
		}
		local formula `"`formula'`token'"'
		gettoken token 0 : 0, parse(",= ") bind
	}
	
	* Preliminary check of the expression by lincom
	cap lincom `formula'
	if _rc!=0 {
		noi di as err "invalid expression"
	}

	* Symbols not allowed in the expression  
	local not_alloed "( ) ^ { } +[ -[ *["
	foreach symb of local not_alloed {
		local check: subinstr local formula  "`symb'" " ", all
		if "`check'"!="`formula'" {
			noi di as err "symbol `symb' not allowed"
			error 198
		}
	}

	* Parse syntax
	local 0 `",`0'"'
	syntax [, Level(cilevel) EForm LONESide RONESide cformat(string)]	

	* Options
	if "`loneside'"!="" & "`roneside'"!="" {
	    noi di as err "Only one of {bf:loneside} and {bf:roneside} can be specified"
		error 198
	}
	if "`loneside'"=="" & "`roneside'"=="" {
		local lev_lb	=(100-`level')/2
		local lev_ub	=100-`lev_lb'
	}
	else {
		local lev_lb	=100-`level'
		local lev_ub	=100-`lev_lb'
	}

	* Estimates from latest wals regression
	tempname b bias V 
	matrix `b'		=e(b)
	matrix `bias'	=e(bias)
	matrix `V'		=e(V)
	local wals_reps `e(reps)'

	* List and number of all covariates 
	local Xlist0 "`e(allvars)'"	
	local k0=e(k0)

	* List and number of non-collinear covariates 
	local Xlist "`e(focvars)' `e(auxvars)'"				
	if "`e(constype)'"!="noconstant" {
		local Xlist: subinstr local Xlist " _cons" ""	
		local Xlist "`Xlist' _cons"	
	}
	local k=e(rank)

	* List and number of omitted covariates 
	local Xlist_omit: list Xlist0 - Xlist
	
	* Tokenize formula 
	local tok_form: subinstr local formula  "+" " +", all	
	local tok_form: subinstr local tok_form "-" " -", all	
	local tok_form: list retokenize tok_form
	
	* Variables involved in the linear combination 
	local vlist: subinstr local tok_form "-" " ", all
	local vlist: subinstr local vlist "+" " ", all
	local vlist: subinstr local vlist "*" " ", all
	local vlist: subinstr local vlist "/" " ", all
	local vlist: subinstr local vlist "_coef[" "", all
	local vlist: subinstr local vlist "_b[" "", all
	local vlist: subinstr local vlist "]" "", all
	local vlist: list retokenize vlist
	local vlist2
	foreach xx of local vlist {
		capture confirm number `xx'
		if _rc!=0 {
			local vlist2 "`vlist2' `xx'"
		}
	}
	local vlist2: list retokenize vlist2
	local vlist: list uniq vlist2
	local check: list vlist in Xlist0
	if `check'==0 {
		foreach vv of local vlist {
			local chech_vv: list vv in Xlist0
			if `chech_vv'==0 {
				noi di as err "[`vv'] not found"
				error 198
			}
		}
		error 198
	}

	* Position of variables in Xlist0
	local pos_vlist 
	foreach vv of local vlist {
		local pos: list posof "`vv'" in Xlist0
		local pos_vlist "`pos_vlist' `pos'"
	}	
	local pos_vlist: list retokenize pos_vlist

	* Vector of coefficients
	tempname C c0
	local c0=0
	matrix `C'=J(1,`k0',0)
	matrix colname `C'= `Xlist0'
	foreach xx of local tok_form {
		* drop `-' and `+'  and determine sign of token
		local yy1: subinstr local xx "-" "", count(local minus)
		local yy1: subinstr local yy1 "+" ""
		if `minus'==1 	local VV_sign "-1"
		else 			local VV_sign "+1"
		
		* transform (`*', `/', _coef[], and _b[]) into " " 
		local yy2: subinstr local yy1 "*" " ", all
		local yy2: subinstr local yy2 "/" " ", all
		local yy2: subinstr local yy2 "_coef[" " ", all
		local yy2: subinstr local yy2 "_b[" " ", all
		local yy2: subinstr local yy2 "]" " ", all
		local yy2: list retokenize yy2

		* drop numbers
		local yy3
		foreach zz of local yy2 {
			capture confirm number `zz'
			if _rc!=0 {
				local yy3 "`yy3' `zz'"
			}
		}
		local yy3: list retokenize yy3
		if "`yy3'"=="" {
		    * put yy2 in the c0 coefficient
		    cap confirm number `yy2'
			if _rc!=0 {
			    noi di as err "invalid expression: check `xx'"
				error 198
			}
			local c0=`c0'`VV_sign'*`yy2'
		}	
		else {
		    * position of the variable in `vlist'
		    local yy3_pos: list posof "`yy3'" in vlist
			
			* put coefficient c_h in the C matrix
			local VV: word `yy3_pos' of `vlist'
			local VV_pos: word `yy3_pos' of `pos_vlist'
			local VV_coef: subinstr local yy1 "_coef[`VV']" "1", count(local isvar)
			if `isvar'==0 {
				local VV_coef: subinstr local yy1 "_b[`VV']" "1", count(local isvar2)
				if `isvar2'==0 {
					local VV_coef: subinstr local yy1 "`VV'" "1"
				}
			}
			matrix `C'[1,`VV_pos']=`C'[1,`VV_pos']`VV_sign'*`VV_coef'
		}
	}
	if "`eform'"!="" & `c0'!=0 {
		noi di as err "additive constant term not allowed with {bf:eform} option"
		exit 198
	}
	tempname C_mata c0_mata
	mata: `c0_mata'=strtoreal(st_local("c0"))
	mata: `C_mata'=st_matrix("`C'")

	* Point estimate and estimaed moments
	tempname lc_est lc_bias lc_V lc_se lc_rmse
	matrix `lc_est'=`c0'+`C'*`b''
	matrix `lc_bias'=`C'*`bias''
	matrix `lc_V'=`C'*`V'*`C''
	matrix `lc_se'=sqrt(`lc_V'[1,1])
	matrix `lc_rmse'=sqrt(`lc_bias'[1,1]^2 + `lc_se'[1,1]^2)

	* List of variables in the dataset of MC replications
	local wals_bc_list ""
	forvalues h=1(1)`k'{
		local wals_bc_list "`wals_bc_list' wals_bc_`h'"
	}

	* Simulation-based confidence intervals of the linear combination
	preserve 
		* Parse to MATA extended dataset of MC replications (including omitted covariates)
		cap use `"`e(bcsimdata)'"', clear
		if _rc!=0 {
			noi di as err "{p 4 4 2}MC replications of bias-corrected wals estimates, "	/*
				*/ "as previously saved by {bf:`e(cmd)'}, do not exist. You need "	/*
				*/ "to refit the model.{p_end}" 
			exit _rc	
		}
		local rwals_bc_list ""
		local oo=1
		local ss=1
		forvalues rr=1(1)`k0'{
			local xx: word `rr' of `Xlist0'
			local xo: word `oo' of `Xlist_omit'
			if "`xx'"!="`xo'" {
				rename wals_bc_`ss' rwals_bc_`rr'
				local ++ss
			}
			else {
				gen double rwals_bc_`rr'=0
				lab var rwals_bc_`rr' "`xo'"
				local ++oo
			}
			local rwals_bc_list "`rwals_bc_list' rwals_bc_`rr'"
		}
		order `rwals_bc_list'
		tempname wals_bc_repdata lc_repdata
		mata: `wals_bc_repdata'=st_data(., "`rwals_bc_list'")
 
		* MC replications of the linear combination
		tempname lc_repdata
		mata: `lc_repdata'=`c0_mata':+`wals_bc_repdata' * `C_mata''
		clear
		qui set obs `wals_reps'
		tempname wals_bc_lc
		qui gen double `wals_bc_lc'=.
		mata: st_store(.,., `lc_repdata')

		* Quantiles  
		tempname lc_ci_lb lc_ci_ub
		_pctile `wals_bc_lc' , percentiles(`lev_lb' `lev_ub')
		scalar `lc_ci_lb'=r(r1)
		scalar `lc_ci_ub'=r(r2)
	restore
	
	* Display linear combination
	noi di _n in gr " (1) ", _continue
	if `c0'!=0 {
	    noi di in ye `c0', _continue
	}
	forvalues hh=1(1)`k0'{
		if `C'[1,`hh']!=0 {
			local VV: word `hh' of `Xlist0'
			local sign=sign(`C'[1,`hh'])
			if `sign'==1 	local sign "+"
			else			local sign "-"
			noi di in ye "`sign' " abs(`C'[1,`hh']) "*`VV'", _continue
		}
	}
	noi di _n "" 
	
	* Check display format 
	if "`cformat'"=="" {
		local cformat "%8.0g"
	}
	else {
		local cfor_len: subinstr local cformat "%" " "
		local cfor_len: subinstr local cfor_len "." " "
		local cfor_len: word 1 of `cfor_len'
		cap confirm integer number `cfor_len'
		if _rc!=0 {
			local cformat "%8.0g"
		}
		else {
			if abs(`cfor_len')>8 {
				local cformat "%8.0g"
			}
		}	
	}

	* Display estimation results
	local COEF_tit="      "
	if "`e(plugin)'"=="ds" {
		local MOM_tit	"  DS   "
		local CI_tit 	"  MC-DS   "
	}
	else {
		local MOM_tit	"  ML   "
		local CI_tit 	"  MC-ML   "
	}
	local LL =10
	local LLp1=`LL'+1
	local p2=`LL'+8	
	local p3=`LL'+18
	local p4=`LL'+28
	local p5=`LL'+38
	local p6=`LL'+53
	if "`eform'"=="" {
		#delimit;
		noi di as text "{hline `LLp1'}{c TT}{hline 63}"													
			_n _col(`LL') "  {c |}"																
			_col(`p2') "`COEF_tit'"                            									
			_col(`p3') "`MOM_tit'"                            									
			_col(`p4') "`MOM_tit'"														 
			_col(`p5') "`MOM_tit'"														 
			_col(`p6') "`CI_tit'" 	                           									
			_n "{ralign `LL':`y'} {c |}" 
			"     Coef.      Bias    Std.Err.    RMSE      [`level'% Conf. Int.]" 
			_n "{hline `LLp1'}{c +}{hline 63}";	
		#delimit cr	
	}
	else {
		#delimit;
		noi di as text "{hline `LLp1'}{c TT}{hline 63}"	
			_n _col(`LL') "  {c |}"																
			_col(`p2') "`COEF_tit'"                            									
			_col(`p3') "`MOM_tit'"                            									
			_col(`p4') "`MOM_tit'"														 
			_col(`p5') "`MOM_tit'"														 
			_col(`p6') "`CI_tit'" 	                           									
			_n "{ralign `LL':`y'} {c |}"															
			"    exp(b)      Bias    Std.Err.    RMSE      [`level'% Conf. Int.]" 
			_n "{hline `LLp1'}{c +}{hline 63}";	
		#delimit cr	
	}
	local p1=`LL'+3	
	local p2=`LL'+16	
	local p3=`LL'+26
	local p4=`LL'+35
	local p5=`LL'+47
	local p6=`LL'+58
	if "`loneside'"=="" & "`roneside'"=="" {
		if "`eform'"=="" {
			#delimit;
			noi di as text "{ralign `LL':(1)} {c |}  "						
					_col(`p1') as res `cformat'  `lc_est'[1,1] 			"  "
					_col(`p2') as res `cformat'  `lc_bias'[1,1]			"  "
					_col(`p3') as res `cformat'  `lc_se'[1,1]			"  "
					_col(`p4') as res `cformat'  `lc_rmse'[1,1]			"  "
					_col(`p5') as res `cformat'  `lc_ci_lb' 			"  "
					_col(`p6') as res `cformat'  `lc_ci_ub';		
			#delimit cr	
		}
		else {
			#delimit;
			noi di as text "{ralign `LL':(1)} {c |}  "						
					_col(`p1') as res `cformat'  exp(`lc_est'[1,1]) 				"  "
					_col(`p2') as res `cformat'  `lc_bias'[1,1]*exp(`lc_est'[1,1])	"  "
					_col(`p3') as res `cformat'  `lc_se'[1,1]*exp(`lc_est'[1,1])	"  "
					_col(`p4') as res `cformat'  `lc_rmse'[1,1]*exp(`lc_est'[1,1])	"  "
					_col(`p5') as res `cformat'  exp(`lc_ci_lb') 			"  "
					_col(`p6') as res `cformat'  exp(`lc_ci_ub');		
			#delimit cr	
		}
	}
	else if "`loneside'"!="" {
		if "`eform'"=="" {
			#delimit;
			noi di as text "{ralign `LL':(1)} {c |}  "						
					_col(`p1') as res `cformat'  `lc_est'[1,1] 			"  "
					_col(`p2') as res `cformat'  `lc_bias'[1,1]			"  "
					_col(`p3') as res `cformat'  `lc_se'[1,1]			"  "
					_col(`p4') as res `cformat'  `lc_rmse'[1,1]			"  "
					_col(`p5') as res `cformat'  . 						"  "
					_col(`p6') as res `cformat'  `lc_ci_ub';		
			#delimit cr	
		}
		else {
			#delimit;
			noi di as text "{ralign `LL':(1)} {c |}  "						
					_col(`p1') as res `cformat'  exp(`lc_est'[1,1]) 				"  "
					_col(`p2') as res `cformat'  `lc_bias'[1,1]*exp(`lc_est'[1,1])	"  "
					_col(`p3') as res `cformat'  `lc_se'[1,1]*exp(`lc_est'[1,1])	"  "
					_col(`p4') as res `cformat'  `lc_rmse'[1,1]*exp(`lc_est'[1,1])	"  "
					_col(`p5') as res `cformat'  . 									"  "
					_col(`p6') as res `cformat'  exp(`lc_ci_ub');		
			#delimit cr	
		}
	}
	else {
		if "`eform'"=="" {
			#delimit;
			noi di as text "{ralign `LL':(1)} {c |}  "						
					_col(`p1') as res `cformat'  `lc_est'[1,1] 			"  "
					_col(`p2') as res `cformat'  `lc_bias'[1,1]			"  "
					_col(`p3') as res `cformat'  `lc_se'[1,1]			"  "
					_col(`p4') as res `cformat'  `lc_rmse'[1,1]			"  "
					_col(`p5') as res `cformat'  `lc_ci_lb' 			"  "
					_col(`p6') as res `cformat'  .;		
			#delimit cr	
		}
		else {
			#delimit;
			noi di as text "{ralign `LL':(1)} {c |}  "						
					_col(`p1') as res `cformat'  exp(`lc_est'[1,1]) 				"  "
					_col(`p2') as res `cformat'  `lc_bias'[1,1]*exp(`lc_est'[1,1])	"  "
					_col(`p3') as res `cformat'  `lc_se'[1,1]*exp(`lc_est'[1,1])	"  "
					_col(`p4') as res `cformat'  `lc_rmse'[1,1]*exp(`lc_est'[1,1])	"  "
					_col(`p5') as res `cformat'  exp(`lc_ci_lb') 			"  "
					_col(`p6') as res `cformat'  .;		
			#delimit cr	
		}
	}
	noi di as text "{hline `LLp1'}{c BT}{hline 63}"		
	
	* Return estimation results
	matrix `C'=`c0', `C'
	local c_list "c0"
	forvalues h=1(1)`k0' {
	    local c_list "`c_list' c`h'"
	}
	matrix coln `C' = `c_list'
	
	return scalar level		=`level'
	if "`loneside'"=="" & "`roneside'"=="" {
		return scalar ub	=`lc_ci_ub'
		return scalar lb 	=`lc_ci_lb'
	}
	else if "`loneside'"!="" {
		return scalar ub	=`lc_ci_ub'
		return scalar lb 	=.
	}
	else if "`loneside'"!="" {
		return scalar ub	=.
		return scalar lb 	=`lc_ci_lb'
	}
	return scalar rmse 		=`lc_rmse'[1,1]
	return scalar se 		=`lc_se'[1,1]
	return scalar bias	 	=`lc_bias'[1,1]
	return scalar estimate 	=`lc_est'[1,1]
	return matrix lc_coef	=`C'
}
end
*-------------------------------------------------------------------------------------------------


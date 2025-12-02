{smcl}
{* *! version 1.0.0  04oct2024}{...}
{vieweralsosee "glmwals postestimation" "help glmwals_postestimation"}{...}
{vieweralsosee "predict" "help glmwals_p"}{...}
{vieweralsosee "lcwals" "help lcwals"}{...}
{vieweralsosee "margwals" "help margwals"}{...}
{vieweralsosee "wals" "help wals"}{...}
{vieweralsosee "hetwals" "help hetwals"}{...}
{vieweralsosee "ar1wals" "help ar1wals"}{...}
{vieweralsosee "xtwals" "help xtwals"}{...}
{vieweralsosee "[LASSO] dslogit" "help dslogit"}{...}
{vieweralsosee "[LASSO] pologit" "help pologit"}{...}
{vieweralsosee "[LASSO] xpologit" "help xpologit"}{...}
{vieweralsosee "[LASSO] dspoisson" "help dspoisson"}{...}
{vieweralsosee "[LASSO] popoisson" "help popoisson"}{...}
{vieweralsosee "[LASSO] xpopoisson" "help xpopoisson"}{...}
{vieweralsosee "[R] glm" "help glm"}{...}
{vieweralsosee "[D] splitsample" "help splitsample"}{...}
{vieweralsosee "[D] vl" "help vl"}{...}
{viewerjumpto "Syntax" "glmwals##syntax"}{...}
{viewerjumpto "Description" "glmwals##description"}{...}
{viewerjumpto "Options" "glmwals##options"}{...}
{viewerjumpto "Examples" "glmwals##examples"}{...}
{viewerjumpto "Stored results" "glmwals##results"}{...}
{viewerjumpto "References" "glmwals##references"}{...}
{p2colset 1 12 14 2}{...}
{p2col:{bf:glmwals} {hline 2}}Weighted-average least squares (WALS) estimation of univariate generalized linear models (GLMs){p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 16 2}{cmd:glmwals} {depvar} [{cmd:(}{help varlist:{it:focvars}}{cmd:)}] {help varlist:{it:auxvars}}
[{it:{help wals##weight:weight}}]
{ifin}
{cmd:[,}
{help wals##optstbl:{it:wals_options}}
{help glmwals##optstbl:{it:glm_options}}]

{p 4 6 2}
where {depvar} is the dependent variable, 
{cmd:(}{help varlist:{it:focvars}}{cmd:)} are the focus regressors, 
{help varlist:{it:auxvars}} are the auxiliary regressors,
{help wals##optstbl:{it:wals_options}} are the basic options of the {help wals:{it:wals}} command,
and 
{help glmwals##optstbl:{it:glm_options}} are the additional options of the {help glmwals:{it:glmwals}} command.
{help varlist:{it:auxvars}} can be omitted only if one specifies the {cmd:auxconstant} option.{p_end}

{marker optstbl}{...}
{synoptset 30 tabbed}{...}
{synopthdr:{it:glm_options}}
{synoptline}
{syntab :GLM setup}
{synopt :{opt f:amily(familyname)}}specify the distribution of {depvar}; default is the {it: Gaussian} distribution{p_end}
{synopt :{opt l:ink(linkname)}}specify the link function; default is the canonical link for the specified {opt family()}{p_end}
{synopt :{opt off:set(ovar)}}include {it:ovar} in the linear predictor with coefficient constrained to 1{p_end}
{synopt :{opt exp:osure(evar)}}include {it:ln(evar)} in the linear predictor with coefficient constrained to 1{p_end}

{syntab :GLM estimation}
{synopt :{opt ltol:erance(#)}}set tolerance for the IRLS estimates; only for {it: re} with i.i.d. errors{p_end}
{synopt :{opt iter0(#)}}set maximum iterations for the IRLS estimates{p_end}
{synopt :{opt ini:t(imethodname)}}set initial values for the mean of {depvar}{p_end}
{synopt :{opt one:step}}specify the one-step WALS estimates{p_end}
{synopt :{opt tol:erance(#)}}set tolerance for the iterative WALS estimates{p_end}
{synopt :{opt iter:ate(#)}}set maximum iterations for the iterative WALS estimates{p_end}

{syntab :Reporting}
{synopt :{opt showi:rls}}show first-step estimates of the unrestricted model{p_end}
{synopt :{opt nolo:g}}suppress iteration log{p_end}
{synopt :{opt eform}}display exponentiated estimates{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{help varlist:{it:focvars}} and {help varlist:{it:auxvars}} may contain factor variables; see {help fvvarlist}.
{p_end}
{p 4 6 2}
{depvar}, {help varlist:{it:focvars}}, and {help varlist:{it:auxvars}} may contain time-series operators; see {help tsvarlist}.
{p_end}
{marker weight}{...}
{p 4 6 2}
{opt fweight}s, {opt iweight}s, and {opt aweight}s are allowed; see {help weight}.
{p_end}
{p 4 6 2}
See {help glmwals_postestimation:glmwals postestimation} for features available after estimation.
{p_end}

{marker description}{...}
{title:Description}

{pstd}
{cmd:glmwals} provides the WALS estimates 
({help glmwals##magnusetal2010:Magnus, Powell, and Pr{c u:}fer 2010})
of univariate generalized linear models (GLMs) such as
gamma, logit, probit, complementary log-log, and Poisson regressions.

{pstd}
Following {help glmwals##delucaetal2018:De Luca, Magnus and Peracchi (2018)},
we linearize the nonlinear system of likelihood equations as in the Newton-Raphson and 
Fisher scoring algorithms.
The main difference with respect to the iteratively reweighted least squares (IRLS) procedure is that,
instead of iterating the optimization until convergence, one performs a WALS regression on the data transformations
obtained in the first iteration.
A single iteration of the estimation procedure leads to the one-step WALS estimator.
To weaken the dependence of this estimator on the starting values, the default estimation method implemented 
by {com: glmwals} is an iterative procedure
that repeatedly updates the starting values using the one-step WALS estimates from the previous iteration
until some convergence criterion is satisfied.

{pstd}
Basic features of the WALS estimation procedure are controlled by the {help wals##optstbl:{it:wals_options}}. 
For additional details see 
{help glmwals##delucaetal2018:De Luca, Magnus and Peracchi (2018)},
and {help glmwals##delucamagnus2024a:De Luca and Magnus (2024a}, {help glmwals##delucamagnus2024b:2024b)}.

{marker options}{...}
{title:glm_options}

{dlgtab:Model setup}

{pstd}
In addition to the model setup options of the {com: wals} command, {com: glmwals} supports the following options:

{phang}
{opt f:amily(familyname)} specifies the distribution of {depvar}.
The available choices for {it: familyname} are

{center: {opt gau:ssian}       {opt nor:mal}         {opt ig:aussian}       {opt in:ormal}}
{center: {opt gam:ma}          {opt b:inomial}       {opt b:ernoulli}       {opt p:oisson}}

{pstd}
where {opt gaussian} is the default distribution,
{opt igaussian} denotes the inverse Gaussian distribution,
{opt normal} is a synonym for {opt gaussian},
and {opt inormal} is a synonym for {opt igaussian}.
The binomial distribution can be specified in three ways:

{center: {opt binomial}       {opt binomial #}       {opt binomial {varname}}}

{pstd}
where {opt binomial 1} is the same as {opt binomial} or {opt bernoulli},
{opt binomial #} defines an integer scalar # for the number of trials,
and
{opt binomial {varname}} allows the number of trials to vary across observations according to
the values of the variable {varname}.
Distribution names are not case sensitive.
By default, the scale parameter for the WALS regressions is set equal to 1 for discrete distributions 
(i.e., {opt binomial} and {opt poisson})
and equal to the Pearson chi-squared statistic divided by the residual degrees of freedom
for continuous distributions (i.e., {opt gaussian}, {opt igaussian}, and {opt gamma}).
To specify alternative values of this parameter one can use the {opt sigma(#)} option (see {help wals##optstbl:{it:wals_options}}).

{phang}
{opt l:ink(linkname)} specifies the link function that maps
the mean outcome into the linear predictor.
The available choices for {it: linkname} are

{center: {opt i:dentity}       {opt log}        {opt l:ogit}            {opt p:robit}        {opt c:loglog}}
{center: {opt logl:og}         {opt logc}       {opt r:eciprocal}       {opt pow:er #}       {opt opo:wer #}}

{pstd} where # is a real scalar.
The default link is the canonical link of each family, that is:

{center:{cmd:family(gaussian)}        {cmd:link(identity)}}
{center:{cmd:family(igaussian)}       {cmd:link(power -2)}}
{center:{cmd:family(binomial)}        {cmd:link(logit)}   }
{center:{cmd:family(poisson)}         {cmd:link(log)}     }
{center:{cmd:family(nbinomial)}       {cmd:link(log)}     }
{center:{cmd:family(gamma)}           {cmd:link(power -1)}}

{pstd} Except for {opt family(nbinomial)} and {opt link(nbinomial)},
{com: glmwals} supports all combinations of {it: familyname} and {it: linkname}
allowed in the {help glm}.
Note, however, that {com: glmwals} aborts with an error when nonstandard combinations
of {it: familyname} and {it: linkname} lead to convergence problems in the preliminary IRLS estimates of the unrestricted model.

{phang}
{opt off:set(ovar)} specifies that the {it: ovar} variable must be included in the linear
predictor with coefficient constrained to be 1.

{phang}
{opt exp:osure(evar)} specifies that the logarithm of the {it: evar} variable must be included in the linear
predictor with coefficient constrained to be 1.

{dlgtab:GLM estimation}

{phang}
{opt ltol:erance(#)} specifies the tolerance for the deviance in the IRLS estimates of the unrestricted model.
The default is {opt ltolerance(1e-6)}.

{phang} 
{opt iter0(#)} specifies the maximum number of iterations for the IRLS estimates of the unrestricted model.
The default is {opt iter0(300)}.

{phang}
{opt ini:t(imethodname)} specifies the initial values for the mean outcome used to initialize the WALS estimation procedure.
We offer four choices for {it: imethodname}:
{it: irls_u} (the default) uses the predicted values obtained from the unrestricted model;
{it: irls_r} uses the predicted values obtained from the fully restricted model;
while {it: irls_0} and {it: irls_0 {varname}} reproduce, respectively,
the default initial values and the {opt mu({varname})} option of the {com: glm} command.

{phang}
{opt one:step} requests the one-step WALS estimates instead of the default iterative WALS estimates.

{phang}
{opt tol:erance(#)} specifies the tolerance for the convergence criterion of the iterative WALS estimates, 
which is the relative difference in the vector of estimated coefficients from one iteration to the next.
The iterative procedure is initialized using the one-step WALS estimates,
and convergence is declared when the criterion is smaller than the specified tolerance level.
The default is {opt tolerance(1e-6)}.

{phang}
{opt iter:ate(#)} specifies the maximum number of iterations for the iterative WALS estimates.
The default is {opt iterate(300)}.
If convergence is not achieved before reaching the specified threshold, then {com: glmwals} displays the current 
results and then aborts with an error.


{dlgtab:Reporting}

{phang}
{opt showi:rls} displays the IRLS estimates of the unrestricted model. 

{phang}
{opt nolo:g} suppresses the display of iteration log in the iterative WALS estimates.
If specified together with the {opt showirls}, this option also suppresses the display of iteration log in the
IRLS estimates of the unrestricted model.

{phang}
{opt eform} displays the exponentiated estimates. 
This option only affects the displayed results and it can be specified during estimation or on replay.
For details see {help nlcom}.

{marker examples}{...}
{title:Examples}

{hline}
{pstd}
Gaussian model{p_end}
{hline}

{pstd}Setup {p_end}
{phang2}{cmd:. use "https://www.stata-press.com/data/r18/breathe", clear}{p_end}
{phang2}{cmd:. run "https://www.stata-press.com/data/r18/no2"}{p_end}
{phang2}{cmd:. display "$cc"}{p_end}
{phang2}{cmd:. display "$fc"}{p_end}
{phang2}{cmd:. rename no2_class NO2}{p_end}

{pstd}WALS estimates{p_end}
{phang2}{cmd:. wals react (NO2) $cc i.($fc), rseed(12345)}{p_end}

{pstd}WALS estimates using the glmwals command{p_end}
{phang2}{cmd:. glmwals react (NO2) $cc i.($fc), rseed(12345)}{p_end}


{hline}
{pstd}
Gamma model with reciprocal link and frequency weights ({help glmwals##hardinhilbe2018:Hardin and Hilbe 2018, Example 6.2}){p_end}
{hline}

{pstd}Setup{p_end}
{phang2}{cmd:. use "http://www.stata-press.com/data/hh4/claims", clear}{p_end}

{pstd}IRLS estimates{p_end}
{phang2}{cmd:. glm y i.pa i.cg i.va i.pa#i.cg i.pa#i.va i.cg#i.va [fw=number], family(gamma) irls}{p_end}

{pstd}WALS estimates{p_end}
{phang2}{cmd:. glmwals y (i.pa i.cg i.va) i.pa#i.cg i.pa#i.va i.cg#i.va [fw=number], family(gamma) prior(laplace) rseed(12345)}{p_end}


{hline}
{pstd}
Binomial models with logit and probit links{p_end}
{hline}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse lbw, clear}{p_end}

{pstd}IRLS estimates of a logit regression{p_end}
{phang2}{cmd:. glm low lwt age i.race smoke ptl ht ui, family(binomial) irls}{p_end}

{pstd}WALS estimates of a logit regression {p_end}
{phang2}{cmd:. glmwals low (lwt) age i.race smoke ptl ht ui, family(binomial) prior(cauchy) rseed(12345)}{p_end}

{pstd}IRLS estimates of a probit regression{p_end}
{phang2}{cmd:. glm low lwt age i.race smoke ptl ht ui, family(binomial) link(probit) irls}{p_end}

{pstd}WALS estimates of a probit regression {p_end}
{phang2}{cmd:. glmwals low (lwt) age i.race smoke ptl ht ui, family(binomial) link(probit) prior(cauchy) rseed(12345)}{p_end}


{hline}
{pstd}
Poisson model{p_end}
{hline}

{pstd}Setup{p_end}
{phang2}{cmd:. use "https://www.stata-press.com/data/r18/breathe", clear}{p_end}
{phang2}{cmd:. do https://www.stata-press.com/data/r18/no2}{p_end}
{phang2}{cmd:. display "$cc"}{p_end}
{phang2}{cmd:. display "$fc"}{p_end}
{phang2}{cmd:. rename no2_class NO2}{p_end}

{pstd}Cross-fit partialing-out lasso estimates of a Poisson regression{p_end}
{phang2}{cmd:. xpopoisson omissions NO2, controls($cc i.($fc)) rseed(12345)}{p_end}

{pstd}WALS estimates of a Poisson regression {p_end}
{phang2}{cmd:. glmwals omissions (NO2) $cc i.($fc), family(poisson) eform noaux rseed(12345)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:glmwals} stores the following in {cmd:e()}: 

{synoptset 23 tabbed}{...}
{p2col 5 23 26 2: Scalars}{p_end}
{synopt :{cmd:e(N)}}number of observations{p_end}
{synopt :{cmd:e(k0)}}number of columns of {cmd:e(b)} (including base levels of factor variables){p_end}
{synopt :{cmd:e(k1)}}number of (not omitted) focus regressors{p_end}
{synopt :{cmd:e(k2)}}number of (not omitted) auxiliary regressors{p_end}
{synopt :{cmd:e(rank)}}rank of {com:e(V)} (also {com:e(k1)+e(k2)}){p_end}
{synopt :{cmd:e(df_r)}}residual degrees of freedom (only if {cmd:e(approach)="fe"}){p_end}
{synopt :{cmd:e(sigma)}}estimated scale parameter{p_end}
{synopt :{cmd:e(quadnpts)}}number of quadrature points for {it:gauss} method{p_end}
{synopt :{cmd:e(quadatol)}}absolute tolerance for {it:adaptive} method{p_end}
{synopt :{cmd:e(quadrtol)}}relative tolerance for {it:adaptive} method{p_end}
{synopt :{cmd:e(reps)}}number of Monte Carlo replications for confidence intervals{p_end}
{synopt :{cmd:e(level)}}confidence level{p_end}
{synopt :{cmd:e(ic)}}number of iterations{p_end}
{synopt :{cmd:e(converged)}}1 if converged, 0 otherwise{p_end}

{p2col 5 23 26 2: Macros}{p_end}
{synopt :{cmd:e(cmdline)}}command as typed{p_end}
{synopt :{cmd:e(cmd)}}{com:glmwals}{p_end}
{synopt :{cmd:e(title)}}title appearing in header{p_end}
{synopt :{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt :{cmd:e(varfunc)}}program used to calculate the variance function{p_end}
{synopt :{cmd:e(varfunct)}}variance title{p_end}
{synopt :{cmd:evarfuncf)}}variance function{p_end}
{synopt :{cmd:e(link)}}program used to calculate the link function{p_end}
{synopt :{cmd:e(linkt)}}link title{p_end}
{synopt :{cmd:e(linkf)}}link function{p_end}
{synopt :{cmd:e(btrials)}}number of binomial trials{p_end}
{synopt :{cmd:e(esttype)}}type of estimates: {it: one-step}, {it: iterative}{p_end}
{synopt :{cmd:e(allvars)}}names of all regressors{p_end}
{synopt :{cmd:e(omitvars)}}names of omitted regressors{p_end}
{synopt :{cmd:e(focvars)}}names of (not omitted) focus regressors{p_end}
{synopt :{cmd:e(auxvars)}}names of (not omitted) auxiliary regressors{p_end}
{synopt :{cmd:e(constype)}}type of constant term: {it:fnoconstant}, {it:focus}, {it:auxiliary}{p_end}
{synopt :{cmd:e(offset)}}linear offset variable{p_end}
{synopt :{cmd:e(wtype)}}weight type (if specified){p_end}
{synopt :{cmd:e(wexp)}}weight expression (if specified){p_end}
{synopt :{cmd:e(prior)}}prior distribution for the population t-ratios{p_end}
{synopt :{cmd:e(quadmethod)}}quadrature method: {it:gauss}, {it:adaptive}, {it:analytic}{p_end}
{synopt :{cmd:e(plugin)}}estimate of sampling moments: {it:ds}, {it:ml}{p_end}
{synopt :{cmd:e(bcsimdata)}}file of MC replications of the bias-corrected WALS estimator{p_end}
{synopt :{cmd:e(properties)}}{it:b V}{p_end}
{synopt :{cmd:e(predict)}}program to implement {com:predict}{p_end}
{synopt :{cmd:e(marginsok)}}predictions allowed by {com:margins} and {com:margwals}{p_end}

{p2col 5 23 26 2: Matrices}{p_end}
{synopt :{cmd:e(b)}}estimated parameter vector{p_end}
{synopt :{cmd:e(V)}}estimated variance matrix of {cmd:e(b)}{p_end}
{synopt :{cmd:e(bias)}}estimated bias vector of {cmd:e(b)}{p_end}
{synopt :{cmd:e(rmse)}}estimated RMSE vector of {cmd:e(b)}{p_end}
{synopt :{cmd:e(MSE)}}estimated MSE matrix of {cmd:e(b)}{p_end}
{synopt :{cmd:e(ci)}}matrix of confidence intervals{p_end}
{synopt :{cmd:e(t_ratios)}}k2-vector of t-ratios in the normal location model{p_end}
{synopt :{cmd:e(wals_pm)}}k2-vector of posterior means in the normal location model{p_end}
{synopt :{cmd:e(wals_wgt)}}k2-vector of WALS weights (diagonal elements of the WALS weight matrix Λ){p_end}
{synopt :{cmd:e(priorpar)}}vector of prior parameters{p_end}

{p2col 5 23 26 2: Functions}{p_end}
{synopt :{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{marker references}{...}
{title:References}

{marker delucamagnus2024a}{...}
{phang}
De Luca, G., and J.R. Magnus. 2024a.
Weighted-Average Least Squares: Improvements and extensions. 
{it:Stata Journal}, forthcoming.
Available at: {browse "https://www.janmagnus.nl/"}

{marker delucamagnus2024b}{...}
{phang}
------. 2024b.
Weighted-Average Least Squares: Beyond the classical linear regression model. 
{it:Stata Journal}, forthcoming.
Available at: {browse "https://www.janmagnus.nl/"}

{marker delucaetal2018}{...}
{phang}
De Luca, G., J.R. Magnus, and F. Peracchi. 2018.
Weighted-average least squares estimation of generalized linear models. 
{it: Journal of Econometrics} 204: 1-17.
{browse "https://doi.org/10.1016/j.jeconom.2017.12.007"}

{marker hardinhilbe2018}{...}
{phang}
Hardin, J. W., and J. M. Hilbe. 2018. 
{it: Generalized Linear Models and Extensions.} 4th ed. 
College Station, TX: Stata Press.

{marker magnusetal2010}{...}
{phang}
Magnus, J.R., O. Powell, and P. Pr{c u:}fer. 2010.
A comparison of two model averaging techniques with an application to growth empirics.
{it:Journal of Econometrics} 154: 139-153.
{browse "https://doi.org/10.1016/j.jeconom.2009.07.004"}
  
{title:Authors}

{pstd}Giuseppe De Luca{p_end}
{pstd}University of Palermo{p_end}
{pstd}Palermo, Italy{p_end}
{pstd}emal:giuseppe.deluca@unipa.it{p_end}

{pstd}Jan R. Magnus {p_end}
{pstd}Vrije Universiteit Amsterdam{p_end}
{pstd}Amsterdam, The Netherlands{p_end}
{pstd}emal:jan@janmagnus.nl{p_end}

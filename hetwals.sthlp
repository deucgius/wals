{smcl}
{* *! version 1.0.0  04oct2024}{...}
{vieweralsosee "hetwals postestimation" "help hetwals_postestimation"}{...}
{vieweralsosee "predict" "help hetwals_p"}{...}
{vieweralsosee "lcwals" "help lcwals"}{...}
{vieweralsosee "margwals" "help margwals"}{...}
{vieweralsosee "wals" "help wals"}{...}
{vieweralsosee "ar1wals" "help ar1wals"}{...}
{vieweralsosee "xtwals" "help xtwals"}{...}
{vieweralsosee "glmwals" "help glmwals"}{...}
{vieweralsosee "[R] hetregress" "help hetregress"}{...}
{vieweralsosee "[D] splitsample" "help splitsample"}{...}
{vieweralsosee "[D] vl" "help vl"}{...}
{viewerjumpto "Syntax" "hetwals##syntax"}{...}
{viewerjumpto "Description" "hetwals##description"}{...}
{viewerjumpto "Options" "hetwals##options"}{...}
{viewerjumpto "Examples" "hetwals##examples"}{...}
{viewerjumpto "Stored results" "hetwals##results"}{...}
{viewerjumpto "References" "hetwals##references"}{...}
{p2colset 1 12 14 2}{...}
{p2col:{bf:hetwals} {hline 2}}Weighted-average least squares (WALS) estimation of linear models with multiplicative heteroskedasticity{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}{cmd:hetwals} {depvar} [{cmd:(}{help varlist:{it:focvars}}{cmd:)}] {help varlist:{it:auxvars}}
[{it:{help wals##weight:weight}}]
{ifin}
{cmd:,}
{cmd:het(}{help varlist:{it:hetvars}}{cmd:)} 
[{help wals##optstbl:{it:wals_options}}
{help hetwals##optstbl:{it:het_options}}]

{p 4 6 2}
where {depvar} is the dependent variable, 
{cmd:(}{help varlist:{it:focvars}}{cmd:)} are the focus regressors, 
{help varlist:{it:auxvars}} are the auxiliary regressors,
{help varlist:{it:hetvars}} are the regressors of the variance function,
{help wals##optstbl:{it:wals_options}} are the basic options of the {help wals:{it:wals}} command,
and 
{help hetwals##optstbl:{it:het_options}} are the additional options of the {help hetwals:{it:hetwals}} command.
{help varlist:{it:auxvars}} can be omitted only if one specifies the {cmd:auxconstant} option.{p_end}

{marker optstbl}{...}
{synoptset 30 tabbed}{...}
{synopthdr:{it:het_options}}
{synoptline}
{syntab :Estimating the variance function}
{synopt :{opt two:step}}use the two-step GLS estimator; default is the ML estimator{p_end}
{synopt :{opt hconst:raints(numlist|matname)}}specifies linear constraints for the variance function{p_end}
{synopt :{opt maximize_options}}control features of the maximization procedure{p_end}

{syntab :Reporting}
{synopt :{opt waldhet}}perform the Wald test for homoskedasticity instead of the likelihood ratio test{p_end}
{synopt :{opt showhetreg}}show first-step estimates of the unrestricted model{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{help varlist:{it:focvars}}, {help varlist:{it:auxvars}}, and {help varlist:{it:hetvars}} may contain factor variables; see {help fvvarlist}.
{p_end}
{p 4 6 2}
{depvar}, {help varlist:{it:focvars}}, {help varlist:{it:auxvars}}, and {help varlist:{it:hetvars}} may contain time-series operators; see {help tsvarlist}.
{p_end}
{marker weight}{...}
{p 4 6 2}
{opt fweight}s, {opt iweight}s, and {opt aweight}s are allowed only when the variance function is estimated by the ML method; see {help weight}.
{p_end}
{p 4 6 2}
The {opt sigma(#)} option from {help wals##optstbl:{it:wals_options}} is not supported.
{p_end}
{p 4 6 2}
See {help hetwals_postestimation:hetwals postestimation} for features available after estimation.
{p_end}

{marker description}{...}
{title:Description}

{pstd}
{cmd:hetwals} computes the WALS estimates of a linear regression model 
with multiplicative heteroskedasticity 
({help hetwals##magnusetal2010:Magnus, Powell, and Pr{c u:}fer 2010}; 
{help hetwals##magnusetal2011:Magnus, Wan, and Zhang 2011}). 
Estimation is based on a feasible
generalized least squares (FGLS) strategy. 
In the first step, we estimate the parameters of the variance function from the unrestricted model 
using  one of the two methods available in the {help hetregress} command:
the (default) ML estimator or Harvey's two-step GLS estimator ({help hetwals##harvey1976:Harvey 1976}).
In the second step, we fit a weighted WALS regression of {depvar} on {help varlist:{it:focvars}} and {help varlist:{it:auxvars}}
with analytic weights equal to the reciprocals of the estimated variances.

{pstd}
Features of the second-step WALS estimates are controlled by the {help wals##optstbl:{it:wals_options}}. 
Also notice that the second-step WALS estimates of the mean function are conditional on the first-step ML/GLS estimates of the variance function.
For additional information on the {help wals} and {help hetwals} commands see {help hetwals##delucamagnus2024a:De Luca and Magnus (2024a,} {help hetwals##delucamagnus2024b:2024b)}.

{marker options}{...}
{title:het_options}

{dlgtab:Estimating the variance function}

{phang} 
{opt two:step} specifies the two-step GLS estimator, instead of the default ML estimator.
This estimation method does not support weights.

{phang}
{opt hconst:raints(numlist|matname)} specifies linear constraints to be applied in
the first-step estimates of the variance function (see {help constraint} and {help makecns}).
This option can only be used with the default ML estimator.

{phang}
{opt maximize_options} include the following options to control the default features of the maximization process (see {help Maximize}): 
{opt tech:nique}({help maximize##algorithm_spec:algorithm_spec}),
{opt tol:erance(#)},
{opt ltol:erance(#)},
{opt nrtol:erance(#)},
{opt nonrtol:erance(#)},
{opt iter:ate(#)},
{opt from}({help maximize##init_specs:init_spec}),
and
{opt dif:ficult}.

{dlgtab:Reporting}

{phang}
{opt waldhet} displays the Wald test for homoskedasticity in the unrestricted model
instead of the likelihood ratio test.

{phang}
{opt showh:etreg} displays the first-step estimates of the unrestricted model. 
For the ML estimator, this option can also be specified together with
the {opt nolog} option to suppress the display of the iteration log.


{marker examples}{...}
{title:Examples}

{hline}
{pstd}Load data{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}

{pstd} ML estimates of a linear model with multiplicative heteroskedasticity {p_end}
{phang2}{cmd:. hetregress price c.mpg##c.mpg length i.foreign c.mpg#c.mpg, het(length i.foreign c.mpg##c.mpg)}{p_end}

{pstd} WALS estimates of a linear model with multiplicative heteroskedasticity {p_end}
{phang2}{cmd:. hetwals price (c.mpg) c.mpg#c.mpg length i.foreign, het(length i.foreign c.mpg##c.mpg) rseed(12345)}{p_end}

{pstd} Heteroskedastic WALS estimates based on the Harvey's two-step GLS estimator of the variance function{p_end}
{phang2}{cmd:. hetwals price (c.mpg) c.mpg#c.mpg length i.foreign, het(length i.foreign c.mpg##c.mpg) twostep showhet rseed(12345)}{p_end}

{pstd} Heteroskedastic WALS estimates based on the horseshoe prior {p_end}
{phang2}{cmd:. hetwals price (c.mpg) c.mpg#c.mpg length i.foreign, het(length i.foreign c.mpg##c.mpg) prior(horseshoe) rseed(12345)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:hetwals} stores the following in {cmd:e()}: 

{synoptset 23 tabbed}{...}
{p2col 5 23 26 2: Scalars}{p_end}
{synopt :{cmd:e(N)}}number of observations{p_end}
{synopt :{cmd:e(k0)}}number of columns of {cmd:e(b)} (including base levels of factor variables){p_end}
{synopt :{cmd:e(k1)}}number of (not omitted) focus regressors (mean function){p_end}
{synopt :{cmd:e(k2)}}number of (not omitted) auxiliary regressors (mean function){p_end}
{synopt :{cmd:e(rank)}}rank of {com:e(V)} (also {com:e(k1)+e(k2)}){p_end}
{synopt :{cmd:e(df_r)}}residual degrees of freedom (mean function){p_end}
{synopt :{cmd:e(quadnpts)}}number of quadrature points for {it:gauss} method{p_end}
{synopt :{cmd:e(quadatol)}}absolute tolerance for {it:adaptive} method{p_end}
{synopt :{cmd:e(quadrtol)}}relative tolerance for {it:adaptive} method{p_end}
{synopt :{cmd:e(reps)}}number of Monte Carlo replications for confidence intervals{p_end}
{synopt :{cmd:e(level)}}confidence level{p_end}
{synopt :{cmd:e(hchi2)}}chi-squared test for homoskedasticity{p_end}
{synopt :{cmd:e(hchi2_df)}}degrees of freedom of {com:e(hchi2)}{p_end}
{synopt :{cmd:e(hchi2_p)}}p-value of {com:e(hchi2)}{p_end}

{p2col 5 23 26 2: Macros}{p_end}
{synopt :{cmd:e(cmdline)}}command as typed{p_end}
{synopt :{cmd:e(cmd)}}{com:hetwals}{p_end}
{synopt :{cmd:e(title)}}title appearing in header{p_end}
{synopt :{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt :{cmd:e(allvars)}}names of all regressors (mean function){p_end}
{synopt :{cmd:e(omitvars)}}names of omitted regressors (mean function){p_end}
{synopt :{cmd:e(focvars)}}names of (not omitted) focus regressors (mean function){p_end}
{synopt :{cmd:e(auxvars)}}names of (not omitted) auxiliary regressors (mean function){p_end}
{synopt :{cmd:e(constype)}}type of constant term in the mean function: {it:noconstant}, {it:focus}, {it:auxiliary}{p_end}
{synopt :{cmd:e(hvars)}}names of regressors in the variance function{p_end}
{synopt :{cmd:e(hmethod)}}estimation method for the variance function: {it:ml}, {it:twostep}{p_end}
{synopt :{cmd:e(htest)}}type of test for \stcmd{e(hchi2)}: LR, Wald{p_end}
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
{synopt :{cmd:e(b)}}estimated parameter vector of the mean function{p_end}
{synopt :{cmd:e(V)}}estimated variance matrix of {cmd:e(b)}{p_end}
{synopt :{cmd:e(bias)}}estimated bias vector of {cmd:e(b)}{p_end}
{synopt :{cmd:e(rmse)}}estimated RMSE vector of {cmd:e(b)}{p_end}
{synopt :{cmd:e(MSE)}}estimated MSE matrix of {cmd:e(b)}{p_end}
{synopt :{cmd:e(ci)}}matrix of confidence intervals (mean function){p_end}
{synopt :{cmd:e(hb)}}estimated parameter vector of the variance function{p_end}
{synopt :{cmd:e(hV)}}estimated variance matrix of {cmd:e(hb)}{p_end}
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

{marker harvey1976}{...}
{phang}
Harvey, A. C. 1976.
Estimating regression models with multiplicative heteroscedasticity. 
{it:Econometrica} 44: 461–465.
{browse "https://doi.org/10.2307/1913974"}

{marker magnusetal2010}{...}
{phang}
Magnus, J.R., O. Powell, and P. Pr{c u:}fer. 2010.
A comparison of two model averaging techniques with an application to growth empirics.
{it:Journal of Econometrics} 154: 139-153.
{browse "https://doi.org/10.1016/j.jeconom.2009.07.004"}
  
{marker magnusetal2011}{...}
{phang}
Magnus, J.R., A.T.K. Wan, and X. Zhang. 2011.
Weighted average least squares estimation with nonspherical disturbances and an application 
to the Hong Kong housing market.
{it:Computational Statistics & Data Analysis} 55: 1331-1341.
{browse "https://doi.org/10.1016/j.csda.2010.09.023"}

{title:Authors}

{pstd}Giuseppe De Luca{p_end}
{pstd}University of Palermo{p_end}
{pstd}Palermo, Italy{p_end}
{pstd}emal:giuseppe.deluca@unipa.it{p_end}

{pstd}Jan R. Magnus {p_end}
{pstd}Vrije Universiteit Amsterdam{p_end}
{pstd}Amsterdam, The Netherlands{p_end}
{pstd}emal:jan@janmagnus.nl{p_end}

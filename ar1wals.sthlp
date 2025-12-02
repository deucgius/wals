{smcl}
{* *! version 1.0.0  04oct2024}{...}
{vieweralsosee "ar1wals postestimation" "help ar1wals_postestimation"}{...}
{vieweralsosee "predict" "help ar1wals_p"}{...}
{vieweralsosee "lcwals" "help lcwals"}{...}
{vieweralsosee "margwals" "help margwals"}{...}
{vieweralsosee "wals" "help wals"}{...}
{vieweralsosee "hetwals" "help hetwals"}{...}
{vieweralsosee "xtwals" "help xtwals"}{...}
{vieweralsosee "glmwals" "help glmwals"}{...}
{vieweralsosee "[TS] prais" "help prais"}{...}
{vieweralsosee "[D] splitsample" "help splitsample"}{...}
{vieweralsosee "[D] vl" "help vl"}{...}
{viewerjumpto "Syntax" "ar1wals##syntax"}{...}
{viewerjumpto "Description" "ar1wals##description"}{...}
{viewerjumpto "Options" "ar1wals##options"}{...}
{viewerjumpto "Examples" "ar1wals##examples"}{...}
{viewerjumpto "Stored results" "ar1wals##results"}{...}
{viewerjumpto "References" "ar1wals##references"}{...}
{p2colset 1 12 14 2}{...}
{p2col:{bf:ar1wals} {hline 2}}Weighted-average least squares (WALS) estimation of linear models with stationary AR(1) errors{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 16 2}{cmd:ar1wals} {depvar} [{cmd:(}{help varlist:{it:focvars}}{cmd:)}] {help varlist:{it:auxvars}}
{ifin}
{cmd:[,}
{help wals##optstbl:{it:wals_options}}
{help ar1wals##optstbl:{it:ar1_options}}]

{p 4 6 2}
where {depvar} is the dependent variable, 
{cmd:(}{help varlist:{it:focvars}}{cmd:)} are the focus regressors, 
{help varlist:{it:auxvars}} are the auxiliary regressors,
{help wals##optstbl:{it:wals_options}} are the basic options of the {help wals:{it:wals}} command,
and 
{help ar1wals##optstbl:{it:ar1_options}} are the additional options of the {help ar1wals:{it:ar1wals}} command.
{help varlist:{it:auxvars}} can be omitted only if one specifies the {cmd:auxconstant} option.{p_end}

{marker optstbl}{...}
{synoptset 30 tabbed}{...}
{synopthdr:{it:ar1_options}}
{synoptline}
{syntab :Estimating the autocorrelation coefficient ρ}
{synopt :{opt rho:type(rhomethod)}}specify the estimation method{p_end}
{synopt :{opt two:step}}stop after the first iteration{p_end}
{synopt :{opt sse:search}}specify the Hildreth-Lu estimation procedure{p_end}
{synopt :{opt optimize_options}}control features of the optimization process{p_end}

{syntab :AR(1) transformations}
{synopt :{opt corc}}use the Cochrane-Orcutt transformation; default is the Prais-Winsten transformation{p_end}

{syntab :Reporting}
{synopt :{opt showar1}}show first-step estimates of the unrestricted model{p_end}
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
Before using {com:ar1wals}, the dataset in memory must be declared as a time-series (see {help tsset}).
{p_end}
{p 4 6 2}
{help varlist:{it:focvars}} and {help varlist:{it:auxvars}} cannot include time-series transformations of {depvar} and weights are not allowed
{p_end}
{p 4 6 2}
See {help ar1wals_postestimation:ar1wals postestimation} for features available after estimation.
{p_end}

{marker description}{...}
{title:Description}

{pstd}
{cmd:ar1wals} computes the WALS estimates of a linear regression model 
with stationary AR(1) errors 
({help ar1wals##magnusetal2010:Magnus, Powell, and Pr{c u:}fer 2010}; 
{help ar1wals##magnusetal2011:Magnus, Wan, and Zhang 2011}). 
Estimation is based on a feasible
generalized least squares (FGLS) strategy. 
In the first step, we estimate the autocorrelation coefficient ρ from the unrestricted model,
using one of the six methods supported by the {help prais} command.
In the second step, we estimate a WALS regression based on the Cochrane-Orcutt or Prais-Winsten transformations of {depvar}, 
{help varlist:{it:focvars}} and {help varlist:{it:auxvars}}.

{pstd}
Features of the second-step WALS estimates are controlled by the {help wals##optstbl:{it:wals_options}}. 
Also notice that the second-step WALS estimates of the mean function are conditional on the first-step estimates of the autocorrelation coefficient ρ.
For additional information on the {help wals} and {help ar1wals} commands see {help ar1wals##delucamagnus2024a:De Luca and Magnus (2024a,} {help ar1wals##delucamagnus2024b:2024b)}.

{marker options}{...}
{title:ar1_options}

{dlgtab:Estimating autocorrelation}

{phang} 
{opt rho:type(rhomethod)} specifies the estimation method for ρ, where the available choices for {it:rhomethod} are:
{synoptset 23 tabbed}{...}{p_end}
{synopt :{opt reg:ress}}LS estimator in a single-lag residual regression{p_end}
{synopt :{opt freg}}LS estimator in a single-lead residual regression{p_end}
{synopt :{opt tsc:orr}}sample autocorrelation of residuals{p_end}
{synopt :{opt dw}}autocorrelation based on Durbin-Watson statistic {p_end}
{synopt :{opt th:eil}}adjusted Theil's sample autocorrelation {p_end}
{synopt :{opt nag:ar}}adjusted Theil-Nagar's autocorrelation {p_end}
{phang} 
The default estimation method is {opt regress}. 
For a description of the various methods see {help prais}. 

{phang}
{opt two:step} specifies that the chosen estimation method for ρ terminates after the first iteration.
By default all methods are iterated until convergence, but the underlying estimators are efficient at each step.

{phang}
{opt sse:search} specifies the Hildreth-Lu estimation procedure, which uses a grid search to find the value of ρ that minimizes 
the sum-of-squared errors in the Cochrane-Orcutt or Prais-Winsten transformation of the unrestricted model.
For additional details see {help prais}. 

{phang}
{opt optimize_options} include two additional options: 
{opt iter:ate(#)} sets the maximum number of iterations in estimating ρ, 
while {opt tol:erance(\num)} sets the tolerance for the convergence of the {com:prais} estimates of the unrestricted model.
The default values of these options are {opt iterate(300)} and {opt tolerance(1e-6)}.

{dlgtab:AR(1) transformations}

{phang}
{opt corc} specifies that the Cochrane-Orcutt transformation must be used instead of the
default Prais-Winsten transformation.
The Prais-Winsten transformation allows preserving the first observation, which is lost in the Cochrane-Orcutt transformation.

{dlgtab:Reporting}

{phang}
{opt showar1} displays the first-step {com: prais} estimates of the unrestricted model. 
The first-step estimate of ρ is always reported in the output's header of {com: ar1wals}.
The {opt showar1} option can also be specified together with {opt nolog} to suppress the display of
the iteration log when estimating the unrestricted model.

{marker examples}{...}
{title:Examples}

{hline}
{pstd}
Prais-Winsten and Cochrane-Orcutt transformations{p_end}
{hline}

{pstd}Setup {p_end}
{phang2}{cmd:. webuse idle, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}

{pstd}Prais-Winsten AR(1) regression {p_end}
{phang2}{cmd:. prais usr idle l.idle}{p_end}

{pstd}WALS Prais-Winsten AR(1) regression {p_end}
{phang2}{cmd:. ar1wals usr (idle) l.idle, rseed(12345)}{p_end}

{pstd}Cochrane-Orcutt AR(1) regression {p_end}
{phang2}{cmd:. prais usr idle l.idle, corc}{p_end}

{pstd}WALS Cochrane-Orcutt AR(1) regression with cauchy prior {p_end}
{phang2}{cmd:. ar1wals usr (idle) l.idle, corc prior(cauchy) rseed(12345)}{p_end}

{hline}
{pstd}
barium data from {help ar1wals##wooldridge2013:Wooldridge (2013, Examples 10.5 and 12.4)}{p_end}
{hline}

{pstd}Setup{p_end}
{phang2}{cmd:. bcuse barium, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}

{pstd}LS estimates{p_end}
{phang2}{cmd:. regress lchnimp lchempi lrtwex lgas befile6 affile6 afdec6}{p_end}

{pstd}Prais-Winsten AR(1) regression{p_end}
{phang2}{cmd:. prais lchnimp lchempi lrtwex lgas befile6 affile6 afdec6}{p_end}

{pstd}WALS Prais-Winsten AR(1) regression with log prior and various estimators of ρ{p_end}
{phang2}{cmd:. ar1wals lchnimp (lchempi lrtwex lgas) befile6 affile6 afdec6, prior(log) rseed(12345)}{p_end}
{phang2}{cmd:. ar1wals lchnimp (lchempi lrtwex lgas) befile6 affile6 afdec6, prior(log) rhot(theil) rseed(12345)}{p_end}
{phang2}{cmd:. ar1wals lchnimp (lchempi lrtwex lgas) befile6 affile6 afdec6, prior(log) sse rseed(12345)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:ar1wals} stores the following in {cmd:e()}: 

{synoptset 23 tabbed}{...}
{p2col 5 23 26 2: Scalars}{p_end}
{synopt :{cmd:e(N)}}number of observations{p_end}
{synopt :{cmd:e(k0)}}number of columns of {cmd:e(b)} (including base levels of factor variables){p_end}
{synopt :{cmd:e(k1)}}number of (not omitted) focus regressors{p_end}
{synopt :{cmd:e(k2)}}number of (not omitted) auxiliary regressors{p_end}
{synopt :{cmd:e(rank)}}rank of {com:e(V)} (also {com:e(k1)+e(k2)}){p_end}
{synopt :{cmd:e(df_r)}}residual degrees of freedom{p_end}
{synopt :{cmd:e(sigma)}}estimated standard deviation of the errors{p_end}
{synopt :{cmd:e(rho)}}first-step estimate of ρ{p_end}
{synopt :{cmd:e(dw_0)}}Durbin-Watson statistic in the (original) unrestricted model{p_end}
{synopt :{cmd:e(dw)}}Durbin-Watson statistic in the (transformed) unrestricted model{p_end}
{synopt :{cmd:e(quadnpts)}}number of quadrature points for {it:gauss} method{p_end}
{synopt :{cmd:e(quadatol)}}absolute tolerance for {it:adaptive} method{p_end}
{synopt :{cmd:e(quadrtol)}}relative tolerance for {it:adaptive} method{p_end}
{synopt :{cmd:e(reps)}}number of Monte Carlo replications for confidence intervals{p_end}
{synopt :{cmd:e(level)}}confidence level{p_end}

{p2col 5 23 26 2: Macros}{p_end}
{synopt :{cmd:e(cmdline)}}command as typed{p_end}
{synopt :{cmd:e(cmd)}}{com:ar1wals}{p_end}
{synopt :{cmd:e(title)}}title appearing in header{p_end}
{synopt :{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt :{cmd:e(allvars)}}names of all regressors{p_end}
{synopt :{cmd:e(omitvars)}}names of omitted regressors{p_end}
{synopt :{cmd:e(focvars)}}names of (not omitted) focus regressors{p_end}
{synopt :{cmd:e(auxvars)}}names of (not omitted) auxiliary regressors{p_end}
{synopt :{cmd:e(constype)}}type of constant term: {it:noconstant}, {it:focus}, {it:auxiliary}{p_end}
{synopt :{cmd:e(rhotype)}}method specified in the {it: rhotype} option{p_end}
{synopt :{cmd:e(rhomethod)}}estimation method for ρ: {it:twostep}, {it:iterated}, {it:SSE search}{p_end}
{synopt :{cmd:e(tranmeth)}}transformation method: {it: corc}, {it: prais}{p_end}
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

{marker wooldridge2013}{...}
{phang}
Wooldridge, J.M. 2013. 
{it: Introductory Econometrics: A Modern Approach.}
5th Edition, South-Western Pub, Mason.

{title:Authors}

{pstd}Giuseppe De Luca{p_end}
{pstd}University of Palermo{p_end}
{pstd}Palermo, Italy{p_end}
{pstd}emal:giuseppe.deluca@unipa.it{p_end}

{pstd}Jan R. Magnus {p_end}
{pstd}Vrije Universiteit Amsterdam{p_end}
{pstd}Amsterdam, The Netherlands{p_end}
{pstd}emal:jan@janmagnus.nl{p_end}

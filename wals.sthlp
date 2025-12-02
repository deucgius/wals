{smcl}
{* *! version 3.0.0  04oct2024}{...}
{vieweralsosee "wals postestimation" "help wals_postestimation"}{...}
{vieweralsosee "predict" "help wals_p"}{...}
{vieweralsosee "lcwals" "help lcwals"}{...}
{vieweralsosee "margwals" "help margwals"}{...}
{vieweralsosee "hetwals" "help hetwals"}{...}
{vieweralsosee "ar1wals" "help ar1wals"}{...}
{vieweralsosee "xtwals" "help xtwals"}{...}
{vieweralsosee "glmwals" "help glmwals"}{...}
{vieweralsosee "[BMA] bmaregress" "help bmaregress"}{...}
{vieweralsosee "[LASSO] dsregress" "help dsregress"}{...}
{vieweralsosee "[LASSO] poregress" "help poregress"}{...}
{vieweralsosee "[LASSO] xporegress" "help xporegress"}{...}
{vieweralsosee "[R] regress" "help regress"}{...}
{vieweralsosee "[D] splitsample" "help splitsample"}{...}
{vieweralsosee "[D] vl" "help vl"}{...}
{vieweralsosee "[R] regress" "help regress"}{...}
{viewerjumpto "Syntax" "wals##syntax"}{...}
{viewerjumpto "Description" "wals##description"}{...}
{viewerjumpto "Options" "wals##options"}{...}
{viewerjumpto "Examples" "wals##examples"}{...}
{viewerjumpto "Stored results" "wals##results"}{...}
{viewerjumpto "References" "wals##references"}{...}
{p2colset 1 12 14 2}{...}
{p2col:{bf:wals} {hline 2}}Weighted-average least squares (WALS) estimation of linear regression models{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}{cmd:wals} {depvar} [{cmd:(}{help varlist:{it:focvars}}{cmd:)}] {help varlist:{it:auxvars}}
{ifin}
[{it:{help wals##weight:weight}}]
[{cmd:,} {help wals##optstbl:{it:options}}]

{p 4 6 2}
where {depvar} is the dependent variable, {cmd:(}{help varlist:{it:focvars}}{cmd:)} is the list of focus regressors, 
and {help varlist:{it:auxvars}} is the list of auxiliary regressors.
{help varlist:{it:auxvars}} can be omitted only if one specifies the {cmd:auxconstant} option.{p_end}

{marker optstbl}{...}
{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{syntab :Model setup}
{synopt :{opt auxcons:tant}}treat the constant term as an auxiliary regressor; default treats the constant term as a focus regressor{p_end}
{synopt :{opt nocons:tant}}suppress the constant term from the model; default includes a constant term as a focus regressor{p_end}
{synopt :{opt sigma(#)}}fix the standard deviation of {depvar} to #; default uses the least squares estimator from the unrestricted model{p_end}

{syntab :Choice of the prior}
{synopt :{opth prior:(wals##priorname:priorname)}}specify the prior distribution for the 'population t-ratio' in the normal location model; default is {cmd:prior({it:pareto})}{p_end}

{syntab :Computing the posterior mean}
{synopt :{opth quadm:ethod(wals##qmetname:qmetname)}}define the integration method for computing the posterior mean; default is {cmd:quadmethod({it:gauss})}{p_end}
{synopt :{opt quadn:pts(#)}}set the number of quadrature points for the {it:gauss} integration method; default is {cmd:quadnpts(500)}{p_end}
{synopt :{opt quade:xt(matname)}}define an external Mata matrix containing the quadrature points for the {it:gauss} integration method{p_end}
{synopt :{opt quadat:ol(#)}}set the absolute tolerance for the {it:adaptive} integration method; default is {cmd:quadatol(1e-9)}{p_end}
{synopt :{opt quadrt:ol(#)}}set the relative tolerance for the {it:adaptive} integration method; default is {cmd:quadrtol(1e-7)}{p_end}

{syntab :Estimating sampling moments (bias and variance)}
{synopt :{opth plug:in(wals##pmetname:pmetname)}}specify the plug-in method for estimating the sampling moments; default is {cmd:plugin({it:ds})}{p_end}
{synopt :{opt patht:ab(path)}}define the directory of the datasets containing the tabulated moments of the posterior mean; default is the current working directory{p_end}
{synopt :{opt fast}}compute point estimates only{p_end}

{syntab :Confidence intervals}
{synopt :{opt lev:el(#)}}set the confidence level; default is {cmd:level(95)}{p_end}
{synopt :{opt rep:s(#)}}set the number of Monte Carlo replications for the confidence intervals; default is {cmd:reps(1000)}{p_end}
{synopt :{opt rseed(#)}}set the random-number seed{p_end}
{synopt :{opt sav:ing}{cmd:(}{it:{help filename:filename}}[{cmd:, replace}]{cmd:)}}save the Monte Carlo replications of the bias-corrected estimator to {it:filename}{cmd:.dta}{p_end}

{syntab :Reporting}
{synopt :{opt nohead:er}}suppress display of output header{p_end}
{synopt :{opt nofoc:us}}suppress display of focus coefficients{p_end}
{synopt :{opt noaux:iliary}}suppress display of auxiliary coefficients{p_end}
{synopt :{opt notab:le}}suppress display of coefficient table{p_end}
{synopt :{opt cformat(%fmt)}}set the display format for statistics in the coefficient table; default is %8.0g{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{help varlist:{it:focvars}} and {help varlist:{it:auxvars}} may contain factor variables; see {help fvvarlist}.
{p_end}
{p 4 6 2}
{depvar}, {help varlist:{it:focvars}}, and {help varlist:{it:auxvars}} may contain time-series operators; see {help tsvarlist}.
{p_end}
{marker weight}{...}
{p 4 6 2}{opt fweight}s, {opt iweight}s, and {opt aweight}s are allowed; see
{help weight}.{p_end}
{p 4 6 2}See {help wals_postestimation:wals postestimation} for features
available after estimation.{p_end}

{synoptset 30}{...}
{marker priorname}{...}
{synopthdr:priorname}
{synoptline}
{synopt :{opt par:eto}}Pareto prior; the default{p_end}
{synopt :{opt lap:lace}}Laplace prior{p_end}
{synopt :{opt sub:botin}}Subbotin prior{p_end}
{synopt :{opt wei:bull}}Weibull prior{p_end}
{synopt :{opt cau:chy}}Cauchy prior{p_end}
{synopt :{opt hor:seshoe}}horseshoe prior{p_end}
{synopt :{opt log}}log prior{p_end}
{synoptline}

{p 4 6 2} Names of the prior distributions are not case sensitive and can be abbreviated as indicated by underlining.

{synoptset 30}{...}
{marker qmetname}{...}
{synopthdr:qmetname}
{synoptline}
{synopt :{opt gauss}}Gauss-Laguerre or Gauss-Legendre integration method; the default{p_end}
{synopt :{opt adaptive}}adaptive integration method (see {help mf_quadrature}){p_end}
{synoptline}

{synoptset 30}{...}
{marker pmetname}{...}
{synopthdr:pmetname}
{synoptline}
{synopt :{opt ds}}plug-in double-shrinkage estimators; the default{p_end}
{synopt :{opt ml}}plug-in maximum likelihood estimators{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:wals} uses the weighted-average least squares estimator proposed by 
{help wals##magnusetal2010:Magnus, Powell, and Pr{c u:}fer (2010)}
to fit a linear regression model with uncertainty about the choice of regressors.
See {help wals##magnusdeluca2016:Magnus and De Luca (2016)} for a review of the WALS approach and 
{help wals##delucaetal2022:De Luca, Magnus and Peracchi (2022,} {help wals##delucaetal2023:2023,} {help wals##delucaetal2024:2024)}
for more recent developments of the WALS theory. 

{pstd}
The statistical framework is the classical linear regression model with two subsets of regressors: 
k1 focus regressors that we want in the model on theoretical or other grounds, 
and k2 auxiliary regressors of which we are less certain.
We assume that k2>0 and that k=k1+k2<n, where n denotes the sample size. 

{pstd}
In this setup, there are 2^k2 possible models that contain all focus regressors and a subset
of the auxiliary regressors.
WALS is a frequentist model averaging estimator that relies on a preliminary semi-orthogonal
transformation of the auxiliary regressors and a Bayesian shrinkage analysis of the associated normal location model.
This model averaging estimator is attractive because it performs well in finite samples,
offers a transparent notion of prior ignorance, and is not restricted to sequences of
nested models. Equally important, it is numerically stable and fast to compute due to
the preliminary semi-orthogonal transformation of the auxiliary regressors.

{pstd}
In its Bayesian shrinkage step, WALS places a prior on the 'population t-ratio' of the 
transformed auxiliary parameters. 
This Bayesian approach is used to obtain an admissible estimator of the location parameter, 
but the quality of this estimator is then assessed in a classical frequentist setup 
(see {help wals##delucaetal2022:De Luca, Magnus and Peracchi 2022,} {help wals##delucaetal2023:2023,} {help wals##delucaetal2024:2024}).

{pstd}
Sampling properties of the posterior mean in the normal location model 
carry over, more or less straightforwardly, to the WALS estimator of the regression parameters. 
In particular, under suitable regurality conditions on the prior, 
the finite-sample distribution of WALS is generally nonnormal and the choice of prior may have sizeable
effects on the estimation bias.
In large samples, this estimator is uniformly {it:sqrt-n} consistent
and its asymptotic distribution is (multivariate) normal only under certain
conditions on the auxiliary parameters.
{help wals##delucaetal2023:De Luca, Magnus and Peracchi (2023)} proposed a simulation method for 
obtaining re-centered and asymmetric WALS confidence intervals based on the bias-corrected posterior mean.

{marker options}{...}
{title:Options}

{dlgtab:Model setup}

{phang} 
{opt auxcons:tant} specifies that the constant term should be treated as an auxiliary regressor rather than as a focus regressor.

{phang}
{opt nocons:tant} removes the constant term from the model. The default procedure includes a constant term as a focus regressor.
Only one of the {com:noconstant} and {com:auxconstant} options can be specified.

{phang}
{opt sigma(#)} fixes the standard deviation {it:sigma} of {depvar} to #, where # is a positive real number.
By default, {it:sigma^2} is estimated by its unbiased least squares estimator from the unrestricted model.
This option is typically used in programs that extend the WALS estimation procedure to more general models.

{dlgtab:Choice of the prior}

{phang}
{opt pri:or(priorname)} specifies the prior distribution for the 'population t-ratio' in the Bayesian shrinkage step of the WALS estimation procedure, 
where {it:priorname} can be one the following: {com:pareto} (the default), 
{com:laplace}, {com:subbotin}, {com:weibull}, {com:cauchy}, {com:horseshoe}, and {com:log}.
Names of the prior distributions are not case sensitive and can be abbreviated as indicated by underlining.
The same prior is automatically assigned to all components of the k2-vector of population t-ratios.
Moreover, the free prior parameters are always fixed to their theoretical values to achieve desirable theoretical properties
such as neutrality, robustness, and minimax regret optimality (see De Luca and Magnus 2024).

{dlgtab:Computing the posterior mean}

{phang}
{opt quadm:ethod(qmetname)} specifies the quadrature integration method for computing the posterior mean under
the {com:subbotin}, {com:weibull}, {com:pareto}, {com:cauchy}, {com:horseshoe}, and {com:log} priors.
The Laplace posterior mean is always computed analytically.
The available options for {it:qmetname} are {com:gauss} (the default) and {com:adaptive},
where {com:gauss} uses the Gauss-Legendre quadrature method for the {com:horseshoe} prior and 
the Gauss-Laguerre quadrature method for the other five priors,
while {com:adaptive} uses the adaptive quadrature method for all priors (see {help mf_quadrature}).
The {com:gauss} method is usually faster than the {com:adaptive} method,
especially when the number of auxiliary regressors or the number of Monte Carlo replications
for the simulation-based confidence intervals are large.

{phang}
{opt quadn:pts(#)} sets the number of quadrature points for the {com:gauss} method to #, where # is a positive integer.
The default is {com:quadnpts(500)}.

{phang}
{opt quade:xt(matname)} specifies the name of a two-column Mata matrix that
contains the quadrature points (first column) and weights (second column) needed for
the {com:gauss} method.
In this case, the number of quadrature points is determined implicitly by the number of rows of {it:matname}.
By default, this matrix is generated internally at any call of {com:wals} based on the {com:gauss} method.
However, in repeated calls to {com:wals} (such as Monte Carlo simulations), 
this option allows reducing substantially the computational burden by generating the
quadrature points and weights externally only once.
The mata library {it:lwalsgaussint.mlib} contains
two functions for generating the external Mata matrix of quadrature points and weights:
{com:gausslegendre(#)} is designed for the horseshoe prior, while
{com:gausslaguerre(#)} is designed for the other five priors,
where # denotes the desired number of quadrature points.
See {help wals##delucamagnus2024:De Luca and Magnus (2024)} for additional details. 

{phang}
{opt quadat:ol(#)} sets the absolute tolerance for the convergence criterion of the {com:adaptive} method to #, 
where # is a positive real number. The default is {com:quadatol(1e-9)}.

{phang}
{opt quadrt:ol(#)} sets the relative tolerance for the convergence criterion of the {com:adaptive} method to #, 
where # is a positive real number. The default is {com:quadrtol(1e-7)}.

{dlgtab:Estimating sampling moments (bias and variance)}

{phang}
{opt plug:in(pmetname)} specifies the plug-in method for estimating the sampling moments of the WALS estimator.
The available choices for {it:pmetname} are {com:ds} for the plug-in double-shrinkage estimators (the default) and {com:ml} 
for the plug-in maximum-likelihood estimators.

{phang}
{opt patht:ab(path)} specifies the directory path where {com:wals} should look for
the datasets containing the tabulations of the sampling bias and variance.
There are seven datasets, named according to the available priors:
{it:wals_pm_tsm_pareto}{cmd:.dta},
{it:wals_pm_tsm_laplace}{cmd:.dta},
{it:wals_pm_tsm_subbotin}{cmd:.dta},
{it:wals_pm_tsm_weibull}{cmd:.dta},
{it:wals_pm_tsm_cauchy}{cmd:.dta},
{it:wals_pm_tsm_horseshoe}{cmd:.dta},
and 
{it:wals_pm_tsm_log}{cmd:.dta}.
By default, {com:wals} assumes these datasets are located in the current working directory
or the ado-path.
If they are located in a different directory, then
this option can be used to specify the appropriate directory path (see {help sysdir} and {help findfile})

{phang}
{opt fast} restricts the calculations of {com:wals} to the point estimates only, excluding sampling moments and confidence intervals.
Specifying {com:fast}, all options related to the sampling moments 
and the confidence intervals become ineffective and
it automatically enables the {com:noheader} and {com:notable} options.

{dlgtab:Confidence intervals}

{phang}
{opt lev:el(#)} specifies the confidence level, where # is a percentage between 10 and 99.99.
The default is {com:level(95)} or the confidence level resulting from {it:c(level)} (see {help level}).
This option can be specified during estimation or on replay.

{phang}
{opt rep:s(#)} sets the number of Monte Carlo replications for the simulation-based confidence intervals
to #, where # is an integer. The default is {com:reps(1000)}.

{phang}
{opt rseed(#)} sets the random-number seed, which is crucial for the reproducibility
of the simulation-based confidence intervals.

{phang}
{opt sav:ing(filename[, replace])}
saves the Monte Carlo replications of the bias-corrected WALS estimator to {it:filename}{cmd:.dta}.
This option can be specified during estimation or on replay.
In both cases, {it:replace} specifies to overwrite {it:filename}{cmd:.dta} if it exists.
When {com:saving} is not specified, {com:wals} saves the simulation results
in a temporary file named {it: STWALS_000001}{cmd:.tmp} for later access by its postestimation commands.
This temporary file is overridden at every run of {com:wals}.
The saved dataset contains one observation for each Monte Carlo replication and a set of variables named 
{it:wals_bc_#} for the not omitted regressors, where # is a progressive index for the regressor number.
The order of the {it:wals_bc_#} variables corresponds to the order in which the regressors are declared in 
{cmd:(}{help varlist:{it:focvars}}{cmd:)} and {help varlist:{it:auxvars}}.
The variable labels of {it:wals_bc_#} contain the names of the corresponding regressors.

{dlgtab:Reporting}

{phang}
{opt nohead:er} suppresses the display of the header at the top of the output.
This option can be specified during estimation or on replay.

{phang}
{opt nofoc:us} suppresses the display of the estimation results for the focus parameters from the coefficient table.
This option can be specified during estimation or on replay.

{phang}
{opt noaux:iliary} suppresses the display of the estimation results for the auxiliary parameters from the coefficient table.
This option can be specified during estimation or on replay.

{phang}
{opt notab:le} suppresses the display of the coefficient table from the output.
This option can be specified during estimation or on replay.

{phang}
{opt cformat(%fmt)} specifies the display format for coefficients, estimated moments, and confidence limits 
in the coefficient table. 
The default format is %8.0g and the maximum format width is 8. 

{marker examples}{...}
{title:Examples}

{hline}
{pstd}Setup
{p_end}
{phang2}{cmd:. use "https://www.stata-press.com/data/r18/breathe", clear}{p_end}
{phang2}{cmd:. run "https://www.stata-press.com/data/r18/no2"}{p_end}
{phang2}{cmd:. display "$cc"}{p_end}
{phang2}{cmd:. display "$fc"}{p_end}
{phang2}{cmd:. rename no2_class NO2}{p_end}

{pstd}WALS estimates based on the default Pareto prior of a linear model for reaction time using the constant term and NO2 as focus regressors
and other control variables as auxiliary regressors{p_end}
{phang2}{cmd:. wals react (NO2) $cc i.($fc), rseed(12345)}{p_end}

{pstd}WALS estimates based on the Weibull prior {p_end}
{phang2}{cmd:. wals react (NO2) $cc i.($fc), rseed(12345) nohead noaux prior(wei)}{p_end}

{pstd}WALS estimates based on the horseshoe prior {p_end}
{phang2}{cmd:. wals react (NO2) $cc i.($fc), rseed(12345) nohead noaux prior(hor)}{p_end}

{pstd}WALS estimates based on the horseshoe prior, computing the posterior mean by the adaptive quadrature method{p_end}
{phang2}{cmd:. wals react (NO2) $cc i.($fc), rseed(12345) nohead noaux prior(hor) quadm(adaptive)}{p_end}

{pstd}WALS estimates based on the horseshoe prior, computing the posterior mean by the Gauss-Legendre method with 1000 quadrature points{p_end}
{phang2}{cmd:. wals react (NO2) $cc i.($fc), rseed(12345) nohead noaux prior(hor) quadnpts(1000)}{p_end}

{pstd}WALS estimates based on the horseshoe prior and plug-in maximum-likelihood estimators of the sampling moments{p_end}
{phang2}{cmd:. wals react (NO2) $cc i.($fc), rseed(12345) nohead noaux prior(hor) plugin(ml)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:wals} stores the following in {cmd:e()}: 

{synoptset 23 tabbed}{...}
{p2col 5 23 26 2: Scalars}{p_end}
{synopt :{cmd:e(N)}}number of observations{p_end}
{synopt :{cmd:e(k0)}}number of columns of {cmd:e(b)} (including base levels of factor variables){p_end}
{synopt :{cmd:e(k1)}}number of (not omitted) focus regressors{p_end}
{synopt :{cmd:e(k2)}}number of (not omitted) auxiliary regressors{p_end}
{synopt :{cmd:e(rank)}}rank of {com:e(V)} (also {com:e(k1)+e(k2)}){p_end}
{synopt :{cmd:e(df_r)}}residual degrees of freedom{p_end}
{synopt :{cmd:e(sigma)}}estimated standard deviation of the errors{p_end}
{synopt :{cmd:e(quadnpts)}}number of quadrature points for {it:gauss} method{p_end}
{synopt :{cmd:e(quadatol)}}absolute tolerance for {it:adaptive} method{p_end}
{synopt :{cmd:e(quadrtol)}}relative tolerance for {it:adaptive} method{p_end}
{synopt :{cmd:e(reps)}}number of Monte Carlo replications for confidence intervals{p_end}
{synopt :{cmd:e(level)}}confidence level{p_end}

{p2col 5 23 26 2: Macros}{p_end}
{synopt :{cmd:e(cmdline)}}command as typed{p_end}
{synopt :{cmd:e(cmd)}}{com:wals}{p_end}
{synopt :{cmd:e(title)}}title appearing in header{p_end}
{synopt :{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt :{cmd:e(allvars)}}names of all regressors{p_end}
{synopt :{cmd:e(omitvars)}}names of omitted regressors{p_end}
{synopt :{cmd:e(focvars)}}names of (not omitted) focus regressors{p_end}
{synopt :{cmd:e(auxvars)}}names of (not omitted) auxiliary regressors{p_end}
{synopt :{cmd:e(constype)}}type of constant term: {it:noconstant}, {it:focus},{it:auxiliary}{p_end}
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
{synopt :{cmd:e(V)}}estimated variance matrix{p_end}
{synopt :{cmd:e(bias)}}estimated bias vector{p_end}
{synopt :{cmd:e(rmse)}}estimated RMSE vector{p_end}
{synopt :{cmd:e(MSE)}}estimated MSE matrix{p_end}
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

{marker delucamagnus2011}{...}
{phang}
De Luca, G., and J.R. Magnus. 2011.
Bayesian model averaging and weighted-average least squares: Equivariance, stability, and numerical issues. 
{it: Stata Journal} 11: 518-544.
{browse "https://journals.sagepub.com/doi/pdf/10.1177/1536867X1201100402"}

{marker delucamagnus2024}{...}
{phang}
------. 2024.
Weighted-Average Least Squares: Improvements and extensions. 
{it:Stata Journal}, forthcoming.
Available at: {browse "https://www.janmagnus.nl/"}

{marker delucaetal2022}{...}
{phang}
De Luca, G., J.R. Magnus, and F. Peracchi. 2022.
Sampling properties of the Bayesian posterior mean with an application to WALS estimation.
{it:Journal of Econometrics} 230: 299–317.
{browse "https://doi.org/10.1016/j.jeconom.2021.04.008"}

{marker delucaetal2023}{...}
{phang}
------. 2023.
Weighted-average least squares (WALS): Confidence and prediction intervals.
{it:Computational Economics} 61: 1637–1664.
{browse "https://link.springer.com/article/10.1007/s10614-022-10255-5"}

{marker delucaetal2024}{...}
{phang}
------. 2024.
Bayesian estimation of the normal location model: A non-standard approach.
{it:Oxford Bulletin of Economics and Statistics}, forthcoming.
Available at: {browse "https://www.janmagnus.nl/"}

{marker magnusdeluca2016}{...}
{phang}
Magnus, J.R., and G. De Luca. 2016.
Weighted-average least squares: A review.
{it:Journal of Economic Surveys} 30: 117-148.
{browse "https://doi.org/10.1111/joes.12094"}

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

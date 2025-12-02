{smcl}
{* *! version 1.0.0  04oct2024}{...}
{vieweralsosee "xtwals postestimation" "help xtwals_postestimation"}{...}
{vieweralsosee "predict" "help xtwals_p"}{...}
{vieweralsosee "lcwals" "help lcwals"}{...}
{vieweralsosee "margwals" "help margwals"}{...}
{vieweralsosee "wals" "help wals"}{...}
{vieweralsosee "hetwals" "help hetwals"}{...}
{vieweralsosee "ar1wals" "help ar1wals"}{...}
{vieweralsosee "glmwals" "help glmwals"}{...}
{vieweralsosee "[XT] xtreg" "help xtreg"}{...}
{vieweralsosee "[XT] xtregar" "help xtregar"}{...}
{vieweralsosee "[D] splitsample" "help splitsample"}{...}
{vieweralsosee "[D] vl" "help vl"}{...}
{viewerjumpto "Syntax" "xtwals##syntax"}{...}
{viewerjumpto "Description" "xtwals##description"}{...}
{viewerjumpto "Options" "xtwals##options"}{...}
{viewerjumpto "Examples" "xtwals##examples"}{...}
{viewerjumpto "Stored results" "xtwals##results"}{...}
{viewerjumpto "References" "xtwals##references"}{...}
{p2colset 1 12 14 2}{...}
{p2col:{bf:xtwals} {hline 2}}Weighted-average least squares (WALS) estimation of fixed-effects and random-effects panel-data models 
with either i.i.d. or stationary AR(1) errors{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 16 2}{cmd:xtwals} {depvar} [{cmd:(}{help varlist:{it:focvars}}{cmd:)}] {help varlist:{it:auxvars}}
[{it:{help wals##weight:weight}}]
{ifin}
{cmd:[,}
{help wals##optstbl:{it:wals_options}}
{help xtwals##optstbl:{it:xt_options}}]

{p 4 6 2}
where {depvar} is the dependent variable, 
{cmd:(}{help varlist:{it:focvars}}{cmd:)} are the focus regressors, 
{help varlist:{it:auxvars}} are the auxiliary regressors,
{help wals##optstbl:{it:wals_options}} are the basic options of the {help wals:{it:wals}} command,
and 
{help xtwals##optstbl:{it:xt_options}} are the additional options of the {help xtwals:{it:xtwals}} command.
{help varlist:{it:auxvars}} can be omitted only if one specifies the {cmd:auxconstant} option.{p_end}

{marker optstbl}{...}
{synoptset 30 tabbed}{...}
{synopthdr:{it:xt_options}}
{synoptline}
{syntab :Model setup}
{synopt :{opt re}}specify the random-effects approach; the default approach is fixed-effects{p_end}
{synopt :{opt ar1}}specify the extended setup with AR(1) errors; the default setup assumes i.i.d. errors{p_end}

{syntab :Estimating variances and autocorrelation}
{synopt :{opt sa}}use the Swamy-Arora estimator of the variance components; only for {it: re} with i.i.d. errors{p_end}
{synopt :{opt rho:type(rhomethod)}}specify the estimation method for the autocorrelation coefficient ρ{p_end}
{synopt :{opt rhof(#)}}set ρ equal to a specific value #{p_end}
{synopt :{opt two:step}}compute the two-step estimate of ρ{p_end}

{syntab :Reporting}
{synopt :{opt showxt:reg}}show first-step estimates of the unrestricted model{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{help varlist:{it:focvars}} and {help varlist:{it:auxvars}} may contain factor variables; see {help fvvarlist}.
{p_end}
{p 4 6 2}
{depvar}, {help varlist:{it:focvars}}, and {help varlist:{it:auxvars}} may contain time-series operators; see {help tsvarlist}.
{p_end}
{p 4 6 2}
Before using {com:xtwals}, the dataset in memory must be declared as a panel dataset (see {help xtset}).
The default setup with i.i.d. error requires declaring at least the panel variable, 
while the extended setup with AR(1) errors requires declaring both the panel and time variables.
{p_end}
{marker weight}{...}
{p 4 6 2}
{opt fweight}s and {opt aweight}s are allowed only in the fixed-effects approach. Furthermore, weights must be constant within units of the panel; see {help weight}.
{p_end}
{p 4 6 2}
{help varlist:{it:focvars}} and {help varlist:{it:auxvars}} cannot include lags of {depvar}.
{p_end}
{p 4 6 2}
See {help xtwals_postestimation:xtwals postestimation} for features available after estimation.
{p_end}

{marker description}{...}
{title:Description}

{pstd}
{cmd:xtwals} provides the WALS estimates 
({help xtwals##magnusetal2010:Magnus, Powell, and Pr{c u:}fer 2010})
of fixed-effects and random-effects panel-data models
with either i.i.d. or AR(1) errors.
The fixed-effects approach is motivated by the fact that the WALS estimator
satisfies the Frisch-Waugh-Lovell theorem when the individual effects are 
treated as (nuisance) focus parameters 
({help xtwals##delucamagnus2024c:De Luca and Magnus 2024c}).
Thus, as with the classical fixed-effects estimator, one can perform
WALS on the familiar within-transformations of the data to wipe out the individual
effects.

{pstd}
The random-effects WALS estimator builds on a special case of the feasible
generalized least squares (FGLS) strategy proposed by 
{help xtwals##magnusetal2011:Magnus, Wan, and Zhang (2011)} 
to perform WALS estimation of linear models with nonspherical errors. 
In this special case, we first estimate the variance components from the unrestricted model, and
then perform a WALS regression based on the FGLS transformations of the original
data to account for the stable equicorrelation exhibited by the one-way errors of the
same unit over time.

{pstd}
In addition to the basic setup with i.i.d. errors, we also use the FGLS transformations 
proposed by 
{help xtwals##bhargavaetal1982:Bhargava, Franzini, and Narendranathan (1982)} 
and
{help xtwals##baltagiwu1999:Baltagi and Wu (1999)} 
to analyze an extended setup of the fixed-effects and random-effects models where 
errors are allowed to follow a stationary AR(1) process 
and observations are allowed to be unequally spaced over time.
For details about the WALS estimation procedure and its extension to panel-data models 
see 
{help xtwals##delucamagnus2024a:De Luca and Magnus (2024a,}
{help xtwals##delucamagnus2024b:2024b,} 
{help xtwals##delucamagnus2024c:2024c)}.

{pstd}
In all cases, basic features of the final WALS estimates are controlled by the {help wals##optstbl:{it:wals_options}}. 
The {cmd:xtwals} command does not support the {opt noconstant} and {opt sigma(#)} options.

{marker options}{...}
{title:xt_options}

{dlgtab:Model setup}

{pstd}
The model setup of {com: xtwals} depends on the following options:

{phang}
{opt re} specifies that the random-effects approach must be used instead of the default fixed-effects approach.

{phang}
{opt ar1} specifies that the extended setup with AR(1) errors must used instead of the default setup with i.i.d.
errors.

{pstd}
The default is the fixed-effects model with i.i.d. errors, which does not require the specification of either 
of the two options.
In contrast, the random-effects model with i.i.d. errors requires the specification of the {opt re} option;
the fixed-effects model with AR(1) errors requires the specification of the {opt ar1} option;
and the random-effects model with AR(1) errors requires the specification of both options.

{pstd}
Each model setup implies certain restrictions on other aspects of the model.
Specifically, in the fixed-effects models,
{help varlist:{it:focvars}} and {help varlist:{it:auxvars}} cannot contain time-invariant regressors,
the {opt auxconstant} option is not allowed,
and the availability of weights is restricted to analytic or frequency weights that are constant within units of the panel.
In the fixed-effects model with AR(1) errors, the only estimation methods for ρ that support weights
are {it: regress} and {it: freg}.
The random-effects models do not support weights.
Finally, in the models with AR(1) errors, {help varlist:{it:focvars}} and {help varlist:{it:auxvars}} 
cannot include the time variable.

{dlgtab:Estimating variances and autocorrelation}

{pstd}
The data transformations used by {com: xtwals} may depend on the
autocorrelation coefficient ρ and the variance components.
These parameters are estimated from the unrestricted model using
the {com: xtreg} command for the basic setups with i.i.d. errors and
the {com: xtregar} command for the extended setups with AR(1) errors.
There are four options related to this aspect of the estimation procedure:

{phang}
{opt sa} uses the small-sample Swamy-Arora estimator of the variance of the random effects
instead of the default consistent estimator.
This option can be specified only in the random-effects model with i.i.d. errors.

{phang} 
{opt rho:type(rhomethod)} specifies the estimation method for ρ, where the available choices for {it:rhomethod} are:
{synoptset 23 tabbed}{...}{p_end}
{synopt :{opt reg:ress}}LS estimator in a single-lag residual regression{p_end}
{synopt :{opt freg}}LS estimator in a single-lead residual regression{p_end}
{synopt :{opt tsc:orr}}sample autocorrelation of residuals{p_end}
{synopt :{opt dw}}autocorrelation based on Durbin-Watson statistic {p_end}
{synopt :{opt th:eil}}adjusted Theil's sample autocorrelation {p_end}
{synopt :{opt nag:ar}}adjusted Theil-Nagar's autocorrelation {p_end}
{synopt :{opt one:step}}Baltagi-Wu's one-step estimator {p_end}
{phang} 
The default estimation method is {opt dw}. 
For details on the various methods see {help xtregar}. 

{phang}
{opt two:step} requests that the first six estimation methods for ρ are stopped after the first iteration.
By default, all methods are iterated until convergence.

{phang}
{opt rhof(#)} requests that ρ is set equal to a value # in the [-1,1] interval.

{pstd}
Note that {opt rhotype}, {opt twostep}, and {opt rhof(#)} must be specified together with {opt ar1};
{opt twostep} cannot be specified together with {opt rhotype(onestep)};
and {opt rhof(\num)} cannot be specified together with {opt rhotype} and {opt twostep}.

{dlgtab:Reporting}

{phang}
{opt showxtreg} displays the preliminary {com:xtreg}/{com:xtregar} estimates of the unrestricted model. 

{marker examples}{...}
{title:Examples}

{hline}
{pstd}
Fixed-effects and random-effects models with i.i.d. errors{p_end}
{hline}

{pstd}Setup {p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. xtset idcode}{p_end}

{pstd}Fixed-effects estimates{p_end}
{phang2}{cmd:. xtreg ln_wage grade c.age c.ttl_exp c.tenure not_smsa i.south i.union i.c_city i.ind_code i.occ_code, fe}{p_end}

{pstd}FE-WALS estimates {p_end}
{phang2}{cmd:. xtwals ln_wage (grade c.age c.ttl_exp c.tenure not_smsa i.south i.union) i.c_city i.ind_code i.occ_code, rseed(12345)}{p_end}

{pstd}Random-effects  estimates {p_end}
{phang2}{cmd:. xtreg ln_wage grade c.age c.ttl_exp c.tenure not_smsa i.south i.union i.c_city i.ind_code i.occ_code, re}{p_end}

{pstd}RE-WALS estimates {p_end}
{phang2}{cmd:. xtwals ln_wage (grade c.age c.ttl_exp c.tenure not_smsa i.south i.union) i.c_city i.ind_code i.occ_code, re rseed(12345)}{p_end}

{hline}
{pstd}
Fixed-effects and random-effects models with AR(1) errors (balanced and equally spaced panel){p_end}
{hline}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse grunfeld, clear}{p_end}
{phang2}{cmd:. xtset company time}{p_end}

{pstd}Fixed-effects FGLS estimates{p_end}
{phang2}{cmd:. xtregar invest c.mvalue##c.mvalue c.kstock##c.kstock c.mvalue#c.kstock, fe}{p_end}

{pstd}FE-WALS estimates{p_end}
{phang2}{cmd:. xtwals invest (c.mvalue c.kstock) c.mvalue#c.mvalue c.kstock#c.kstock c.mvalue#c.kstock, ar1 prior(subbotin) rseed(12345)}{p_end}

{pstd}Random-effects FGLS estimates{p_end}
{phang2}{cmd:. xtregar invest c.mvalue##c.mvalue c.kstock##c.kstock c.mvalue#c.kstock, re}{p_end}

{pstd}RE-WALS estimates {p_end}
{phang2}{cmd:. xtwals invest (c.mvalue c.kstock) c.mvalue#c.mvalue c.kstock#c.kstock c.mvalue#c.kstock, re ar1 prior(subbotin) rseed(12345)}{p_end}


{hline}
{pstd}
Fixed-effects and random-effects models with AR(1) errors (unbalanced and unequally spaced panel){p_end}
{hline}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. xtset idcode year}{p_end}
{phang2}{cmd:. xtdes}{p_end}

{pstd}FE-WALS estimates{p_end}
{phang2}{cmd:. xtwals ln_wage (grade c.age c.ttl_exp c.tenure not_smsa i.south i.union) i.c_city i.ind_code i.occ_code, ar1 rhotype(onestep) rseed(12345)}{p_end}

{pstd}RE-WALS estimates {p_end}
{phang2}{cmd:. xtwals ln_wage (grade c.age c.ttl_exp c.tenure not_smsa i.south i.union) i.c_city i.ind_code i.occ_code, re ar1 rhotype(onestep) rseed(12345)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:xtwals} stores the following in {cmd:e()}: 

{synoptset 23 tabbed}{...}
{p2col 5 23 26 2: Scalars}{p_end}
{synopt :{cmd:e(N)}}number of observations{p_end}
{synopt :{cmd:e(n)}}number of units{p_end}
{synopt :{cmd:e(Tmin)}}minimum number of observations across units{p_end}
{synopt :{cmd:e(Tavg)}}average number of observations across units{p_end}
{synopt :{cmd:e(Tmax)}}maximum number of observations across units{p_end}
{synopt :{cmd:e(Tcon)}}1 if observations are constant across units{p_end}
{synopt :{cmd:e(k0)}}number of columns of {cmd:e(b)} (including base levels of factor variables){p_end}
{synopt :{cmd:e(k1)}}number of (not omitted) focus regressors{p_end}
{synopt :{cmd:e(k2)}}number of (not omitted) auxiliary regressors{p_end}
{synopt :{cmd:e(rank)}}rank of {com:e(V)} (also {com:e(k1)+e(k2)}){p_end}
{synopt :{cmd:e(df_r)}}residual degrees of freedom (only if {cmd:e(approach)="fe"}){p_end}
{synopt :{cmd:e(sigma)}}estimated standard deviation of the errors in the WALS regression on the transformed data{p_end}
{synopt :{cmd:e(sig_nu)}}estimated standard deviation of random effects (only if {cmd:e(approach)="re"}){p_end}
{synopt :{cmd:e(sig_e)}}estimated standard deviation of idiosyncratic regression errors (only if {cmd:e(approach)="re"}){p_end}
{synopt :{cmd:e(alpha)}}first-step estimate of the GLS coefficient α (only if {cmd:e(approach)="re"} and {cmd:e(Tcon)=1}){p_end}
{synopt :{cmd:e(rho)}}first-step estimate of the autocorrelation coefficient ρ{p_end}
{synopt :{cmd:e(quadnpts)}}number of quadrature points for {it:gauss} method{p_end}
{synopt :{cmd:e(quadatol)}}absolute tolerance for {it:adaptive} method{p_end}
{synopt :{cmd:e(quadrtol)}}relative tolerance for {it:adaptive} method{p_end}
{synopt :{cmd:e(reps)}}number of Monte Carlo replications for confidence intervals{p_end}
{synopt :{cmd:e(level)}}confidence level{p_end}

{p2col 5 23 26 2: Macros}{p_end}
{synopt :{cmd:e(cmdline)}}command as typed{p_end}
{synopt :{cmd:e(cmd)}}{com:xtwals}{p_end}
{synopt :{cmd:e(title)}}title appearing in header{p_end}
{synopt :{cmd:e(model)}}model setup: {it:fe-iid}, {it:re-iid}, {it:fe-ar1}, {it:re-ar1}{p_end}
{synopt :{cmd:e(approach)}}approach: {it:fe}, {it:re}{p_end}
{synopt :{cmd:e(errors)}}regression errors: {it:iid}, {it:ar1}{p_end}
{synopt :{cmd:e(ivar)}}variable denoting units{p_end}
{synopt :{cmd:e(tvar)}}variable denoting time within units (only if {com: e(errors)="ar1"}){p_end}
{synopt :{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt :{cmd:e(allvars)}}names of all regressors{p_end}
{synopt :{cmd:e(omitvars)}}names of omitted regressors{p_end}
{synopt :{cmd:e(focvars)}}names of (not omitted) focus regressors{p_end}
{synopt :{cmd:e(auxvars)}}names of (not omitted) auxiliary regressors{p_end}
{synopt :{cmd:e(constype)}}type of constant term: {it:focus}, {it:auxiliary}{p_end}
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
{synopt :{cmd:e(alpha_stat)}}summary statistics on GLS coefficients α_i (only if {cmd:e(approach)="re"} and {cmd:e(Tcon)=0}){p_end}
{synopt :{cmd:e(t_ratios)}}k2-vector of t-ratios in the normal location model{p_end}
{synopt :{cmd:e(wals_pm)}}k2-vector of posterior means in the normal location model{p_end}
{synopt :{cmd:e(wals_wgt)}}k2-vector of WALS weights (diagonal elements of the WALS weight matrix Λ){p_end}
{synopt :{cmd:e(priorpar)}}vector of prior parameters{p_end}

{p2col 5 23 26 2: Functions}{p_end}
{synopt :{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{marker references}{...}
{title:References}

{marker baltagiwu1999}{...}
{phang}
Baltagi, B. H., and P. X. Wu. 1999. 
Unequally spaced panel data regressions with AR(1) disturbances. 
{it: Econometric Theory} 15: 814–823.
{browse "https://doi.org/10.1017/S0266466699156020 "}

{marker bhargavaetal1982}{...}
{phang}
Bhargava, A., L. Franzini, and W. Narendranathan. 1982. 
Serial correlation and the fixed effects model. 
{it: Review of Economic Studies} 49: 533–549.
{browse "https://doi.org/10.2307/2297285"}

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

{marker delucamagnus2024c}{...}
{phang}
------. 2024c.
Weighted-average least squares estimation of panel-data models. 
In 
{it: Advances in Shrinkage and Penalized Estimation Strategies: Honoring the Contributions
of A. K. Md. Ehsanes Saleh}, 
ed. M. Arashi and M. Norouzirad, 00–00. New York:
Springer Series "Emerging Topics in Statistics and Biostatistics".
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

{title:Authors}

{pstd}Giuseppe De Luca{p_end}
{pstd}University of Palermo{p_end}
{pstd}Palermo, Italy{p_end}
{pstd}emal:giuseppe.deluca@unipa.it{p_end}

{pstd}Jan R. Magnus {p_end}
{pstd}Vrije Universiteit Amsterdam{p_end}
{pstd}Amsterdam, The Netherlands{p_end}
{pstd}emal:jan@janmagnus.nl{p_end}

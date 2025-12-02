{smcl}
{* *! version 1.0.0  04oct2024}{...}
{vieweralsosee "xtwals" "help xtwals"}{...}
{vieweralsosee "xtwals postestimation" "help xtwals_postestimation"}{...}
{vieweralsosee "lcwals" "help lcwals"}{...}
{vieweralsosee "margwals" "help margwals"}{...}
{vieweralsosee "wals" "help wals"}{...}
{vieweralsosee "hetwals" "help hetwals"}{...}
{vieweralsosee "ar1wals" "help ar1wals"}{...}
{vieweralsosee "glmwals" "help glmwals"}{...}
{viewerjumpto "Syntax" "xtwals_p##syntax"}{...}
{viewerjumpto "Description" "xtwals_p##description"}{...}
{viewerjumpto "Examples" "xtwals_p##examples"}{...}
{viewerjumpto "References" "xtwals_p##references"}{...}
{p2colset 1 12 14 2}{...}
{p2col:{bf:predict} {hline 2}}Weighted-average least squares (WALS) predictions for fixed-effects and random-effects panel-data models 
with either i.i.d. or stationary AR(1) errors{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Syntax of the {com:predict} command for {helpb xtwals} is presented under the following headings:

{phang2}{help xtwals_p##syntax1:Compute linear predictions, moments of linear predictions, and residuals}{p_end}
{phang2}{help xtwals_p##syntax2:Compute confidence and prediction intervals}{p_end}


{marker syntax1}{...}
{phang}{ul:{bf:Compute linear predictions, moments of linear predictions, and residuals}}

{p 8 16 2}
{cmd:predict} {dtype} {newvar} {ifin} 
[{cmd:,} {help xtwals_p##opts1:{it:predict_options1}}]

{marker syntax2}{...}
{phang}{ul:{bf:Compute confidence and prediction intervals}}

{p 8 16 2}
{cmd:predict} {dtype} {newvar}_{it:l}  {newvar}_{it:u} {ifin} 
[{cmd:,} {help xtwals_p##opts2:{it:predict_options2}}]


{synoptset 20 tabbed}{...}
{marker opts1}{...}
{synopthdr :predict_options1}
{synoptline}
{synopt :{opt xb}}linear prediction (default){p_end}
{synopt :{opt biasp}}bias of the linear prediction{p_end}
{synopt :{opt stdp}}standard error of the linear prediction{p_end}
{synopt :{opt rmsep}}root mean squared error of the linear prediction{p_end}
{synopt :{opt bcp}}bias-corrected linear prediction{p_end}
{synopt :{opt stdbcp}}standard error of the bias-corrected linear prediction{p_end}
{synopt :{opt ie}}individual effects{p_end}
{synopt :{opt bcie}}bias-corrected individual effects{p_end}
{synopt :{opt r:esiduals}}residuals from the linear prediction{p_end}
{synopt :{opt bcr:esiduals}}residuals from the bias-corrected linear prediction{p_end}
{synopt :{opt wp}}{it:s} periods-ahead WALS predictions{p_end}
{synopt :{opt bcwp}}{it:s} periods-ahead bias-corrected WALS predictions{p_end}
{synoptline}
{p2colreset}{...}


{synoptset 20 tabbed}{...}
{marker opts2}{...}
{synopthdr :predict_options2}
{synoptline}
{synopt :{opt cint:erval}}confidence intervals for E(y_{i,T+s}|x_{i,T+s}) {p_end}
{synopt :{opt pint:erval}}prediction intervals for y_{i,T+s}|x_{i,T+s}{p_end}
{synopt :{opt lev:el(#)}}sets the confidence level{p_end}
{synopt :{opt rseed(#)}}sets the random-number seed for the prediction intervals{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:predict} uses the current estimation results produced by {helpb xtwals} to compute weighted-average least squares (WALS) predictions 
such as linear predictions, bias-corrected predictions, their sampling moments (bias, SE, and RMSE), individual effects, residuals, 
bias-corrected residuals, {it:s} periods-ahead WALS predictions, confidence intervals for E(y_{i,T+s}|x_{i,T+s}), 
and prediction intervals for y_{i,T+s}|x_{i,T+s}. 

{pstd}
In particular, the first syntax creates in {newvar} one of the statistics available in {help xtwals_p##opts1:{it:predict_options1}}.
If no option is specified, the default {com: xb} option is assumed. 
Estimated sampling moments and bias-corrected statistics depend on the plug-in estimators specified in the last {helpb xtwals} regression 
(see {help wals##options:{it:wals_options}}).
The {opt ie}, {opt bcie}, {opt r:esiduals}, {opt bcr:esiduals}, {opt wp}, and {opt bcwp} options are currently restricted 
to the fixed-effects and random-effects models with i.i.d. errors.
The only statistic allowed with the {help margwals} command is the default linear prediction.

{pstd}
The second syntax creates in {newvar}_{it:l} and {newvar}_{it:u} 
the lower and upper bounds of the confidence intervals for E(y_{i,T+s}|x_{i,T+s}) if one specifies the {cmd:cinterval} option,
and the lower and upper bounds of the prediction intervals for y_{i,T+s}|x_{i,T+s} if one specifies the {cmd:pinterval} option.
Confidence and prediction intervals are currently resticted to the fixed-effects and random-effects models with i.i.d. errors.
The confidence level of both intervals is controlled by the {opt level(#)} option.

{pstd}
Except for the {opt xb}, {opt ie}, and {opt wp} options, 
all statistics require that the last {helpb xtwals} regression was fitted without the {opt fast} option (see {help wals##options:{it:wals_options}}).
By default all statistics are calculated out of sample; type {cmd:predict ... if e(sample) ...} 
if calculations must be restricted to the estimation sample.
For additional details see {help xtwals_p##delucamagnus2024a:De Luca and Magnus (2024a,} {help xtwals_p##delucamagnus2024b:2024b)}.

{marker options}{...}
{title:Options}

{dlgtab:predict_options1}

{phang} 
{opt xb}, the default, calculates the linear prediction without including the individual effects.

{phang} 
{opt biasp} calculates the plug-in estimate of the bias of the linear prediction.

{phang} 
{opt stdp} calculates the plug-in estimate of the standard error of the linear prediction.

{phang} 
{opt rmsep} calculates the plug-in estimate of the root mean squared error of the linear prediction.

{phang} 
{opt bcp} calculates the bias-corrected linear prediction.

{phang} 
{opt stdbcp} calculates the standard error of the bias-corrected linear prediction.

{phang} 
{opt ie} calculates the estimates of the individual fixed-effects or random-effects.
 
{phang} 
{opt bcie} calculates the bias-corrected WALS estimates of the individual effects.

{phang} 
{opt r:esiduals} calculates the residuals from the linear prediction and the individual effects.

{phang} 
{opt bcr:esiduals} calculates the residuals from the bias-corrected linear prediction and the bias-corrected individual effects.

{phang} 
{opt wp} calculates the {it:s} periods-ahead WALS predictions.

{phang} 
{opt bcwp} calculates the {it:s} periods-ahead bias-corrected WALS predictions.

{dlgtab:predict_options2}

{phang} 
{opt cint:erval} calculates the confidence intervals for E(y_{i,T+s}|x_{i,T+s}).

{phang} 
{opt pint:erval} calculates the prediction intervals for y_{i,T+s}|x_{i,T+s}. 

{phang} 
{opt lev:el(#)} specifies the confidence level of the confidence/prediction intervals, where # is a percentage between 10 and 99.99.
The default is {com:level(95)} or the confidence level resulting from {it:c(level)} (see {helpb level}).

{phang} 
{opt rseed(#)} sets the random-number seed to ensure reproducibility of the prediction intervals. 

{marker examples}{...}
{title:Examples}

{hline}
{pstd}
Fixed-effects and random-effects models with i.i.d. errors {p_end}
{hline}

{pstd}
Setup
{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. xtset idcode}{p_end}

{pstd}
Fixed-effects estimates of individual effects
{p_end}
{phang2}{cmd:. xtreg ln_wage grade c.age c.ttl_exp c.tenure not_smsa i.south i.union i.c_city i.ind_code i.occ_code, fe}{p_end}
{phang2}{cmd:. predict xtreg_fe, u}{p_end}

{pstd}
FE-WALS estimates of individual effects
{p_end}
{phang2}{cmd:. xtwals ln_wage (grade c.age c.ttl_exp c.tenure not_smsa i.south i.union) i.c_city i.ind_code i.occ_code, rseed(12345)}{p_end}
{phang2}{cmd:. predict xtwals_fe, ie}{p_end}

{pstd}
Random-effects estimates of individual effects
{p_end}
{phang2}{cmd:. xtreg ln_wage grade c.age c.ttl_exp c.tenure not_smsa i.south i.union i.c_city i.ind_code i.occ_code, re}{p_end}
{phang2}{cmd:. predict xtreg_re, u}{p_end}

{pstd}
RE-WALS estimates of individual effects
{p_end}
{phang2}{cmd:. xtwals ln_wage (grade c.age c.ttl_exp c.tenure not_smsa i.south i.union) i.c_city i.ind_code i.occ_code, re rseed(12345)}{p_end}
{phang2}{cmd:. predict xtwals_re, ie}{p_end}


{hline}
{pstd}
Fixed-effects and random-effects models with AR(1) errors (balanced and equally spaced panel) {p_end}
{hline}

{pstd}
Setup
{p_end}
{phang2}{cmd:. webuse grunfeld, clear}{p_end}
{phang2}{cmd:. xtset company time}{p_end}
 
{pstd}
Linear predictions from FE-WALS estimates
{p_end}
{phang2}{cmd:. xtwals invest (c.mvalue c.kstock) c.mvalue#c.mvalue c.kstock#c.kstock c.mvalue#c.kstock, ar1 prior(subbotin) rseed(12345)}{p_end}
{phang2}{cmd:. predict double xtwals_fe_lp}{p_end}

{pstd}
Linear predictions from  RE-WALS estimates
{p_end}
{phang2}{cmd:. xtwals invest (c.mvalue c.kstock) c.mvalue#c.mvalue c.kstock#c.kstock c.mvalue#c.kstock, re ar1 prior(subbotin) rseed(12345)}{p_end}
{phang2}{cmd:. predict double xtwals_re_lp}{p_end}


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

{title:Authors}

{pstd}Giuseppe De Luca{p_end}
{pstd}University of Palermo{p_end}
{pstd}Palermo, Italy{p_end}
{pstd}emal:giuseppe.deluca@unipa.it{p_end}

{pstd}Jan R. Magnus {p_end}
{pstd}Vrije Universiteit Amsterdam{p_end}
{pstd}Amsterdam, The Netherlands{p_end}
{pstd}emal:jan@janmagnus.nl{p_end}

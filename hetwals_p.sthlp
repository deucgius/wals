{smcl}
{* *! version 1.0.0  04oct2024}{...}
{vieweralsosee "hetwals" "help hetwals"}{...}
{vieweralsosee "hetwals postestimation" "help hetwals_postestimation"}{...}
{vieweralsosee "lcwals" "help lcwals"}{...}
{vieweralsosee "margwals" "help margwals"}{...}
{vieweralsosee "wals" "help wals"}{...}
{vieweralsosee "ar1wals" "help ar1wals"}{...}
{vieweralsosee "xtwals" "help xtwals"}{...}
{vieweralsosee "glmwals" "help glmwals"}{...}
{viewerjumpto "Syntax" "hetwals_p##syntax"}{...}
{viewerjumpto "Description" "hetwals_p##description"}{...}
{viewerjumpto "Examples" "hetwals_p##examples"}{...}
{viewerjumpto "References" "hetwals_p##references"}{...}
{p2colset 1 12 14 2}{...}
{p2col:{bf:predict} {hline 2}}Weighted-average least squares (WALS) predictions for linear models with multiplicative heteroskedasticity{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Syntax of the {com:predict} command for {helpb hetwals} is presented under the following headings:

{phang2}{help hetwals_p##syntax1:Compute linear predictions, moments of linear predictions, and residuals}{p_end}
{phang2}{help hetwals_p##syntax2:Compute confidence and prediction intervals}{p_end}


{marker syntax1}{...}
{phang}{ul:{bf:Compute linear predictions, moments of linear predictions, and residuals}}

{p 8 16 2}
{cmd:predict} {dtype} {newvar} {ifin} 
[{cmd:,} {help hetwals_p##opts1:{it:predict_options1}}]

{marker syntax2}{...}
{phang}{ul:{bf:Compute confidence and prediction intervals}}

{p 8 16 2}
{cmd:predict} {dtype} {newvar}_{it:l}  {newvar}_{it:u} {ifin} 
[{cmd:,} {help hetwals_p##opts2:{it:predict_options2}}]


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
{synopt :{opt r:esiduals}}residuals from the linear prediction{p_end}
{synopt :{opt bcr:esiduals}}residuals from the bias-corrected linear prediction{p_end}
{synopt :{opt sigma}}estimated standard deviations of the errors{p_end}
{synoptline}
{p2colreset}{...}


{synoptset 20 tabbed}{...}
{marker opts2}{...}
{synopthdr :predict_options2}
{synoptline}
{synopt :{opt cint:erval}}confidence intervals for E(y_i|x_i) {p_end}
{synopt :{opt pint:erval}}prediction intervals for y_i|x_i{p_end}
{synopt :{opt lev:el(#)}}sets the confidence level{p_end}
{synopt :{opt rseed(#)}}sets the random-number seed for the prediction intervals{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:predict} uses the current estimation results produced by {helpb hetwals} to compute weighted-average least squares (WALS) predictions 
such as linear predictions, bias-corrected predictions, their sampling moments (bias, SE, and RMSE), residuals, bias-corrected residuals,
estimated standard deviations of the errors, confidence intervals for E(y_i|x_i), and prediction intervals for y_i|(x_i,v_i). 

{pstd}
In particular, the first syntax creates in {newvar} one of the statistics available in {help hetwals_p##opts1:{it:predict_options1}}.
If no option is specified, the default {com: xb} option is assumed. 
Estimated sampling moments and bias-corrected statistics depend on the plug-in estimators specified in the last {helpb hetwals} regression 
(see {help wals##options:{it:wals_options}}). 
The first eight statistics 
(i.e., {opt xb}, {opt biasp}, {opt stdp}, {opt rmsep}, {opt bcp}, {opt stdbcp}, {opt r:esiduals}, and {opt bcr:esiduals}) 
are conditional on the first-step estimates of the variance function.
The only statistic allowed with the {help margwals} command is the default linear prediction.

{pstd}
The second syntax creates in {newvar}_{it:l} and {newvar}_{it:u} 
the lower and upper bounds of the confidence intervals for E(y_i|x_i) if one specifies the {cmd:cinterval} option,
and the lower and upper bounds of the prediction intervals for y_i|(x_i,v_i) if one specifies the {cmd:pinterval} option.
The confidence level of both intervals is controlled by the {opt level(#)} option.

{pstd}
Except for the linear prediction, all statistics require that 
the last {helpb hetwals} regression was fitted without the {opt fast} option (see {help wals##options:{it:wals_options}}).
All statistics are available both in sample and out of sample (default); type {cmd:predict ... if e(sample) ...} 
if calculations must be restricted to the estimation sample.
For additional details see {help hetwals##delucamagnus2024a:De Luca and Magnus (2024a,} {help hetwals##delucamagnus2024b:2024b)}.

{marker options}{...}
{title:Options}

{dlgtab:predict_options1}

{phang} 
{opt xb}, the default, calculates the linear prediction (i.e., the WALS predictor).

{phang} 
{opt biasp} calculates the bias of the linear prediction.

{phang} 
{opt stdp} calculates the standard error of the linear prediction.

{phang} 
{opt rmsep} calculates the root mean squared error of the linear prediction.

{phang} 
{opt bcp} calculates the bias-corrected linear prediction.

{phang} 
{opt stdbcp} calculates the standard error of the bias-corrected linear prediction.

{phang} 
{opt r:esiduals} calculates the residuals from the linear prediction.

{phang} 
{opt bcr:esiduals} calculates the residuals from the bias-corrected linear prediction.

{phang} 
{opt sigma} calculates the estimated standard deviations of the errors.

{dlgtab:predict_options2}

{phang} 
{opt cint:erval} calculates the confidence intervals for E(y_i|x_i).

{phang} 
{opt pint:erval} calculates the prediction intervals for y_i|(x_i,v_i), where x_i are the (focus and auviliary) regressors of the mean function 
and v_i are the regressors of the variance function. 

{phang} 
{opt lev:el(#)} specifies the confidence level of the confidence/prediction intervals, where # is a percentage between 10 and 99.99.
The default is {com:level(95)} or the confidence level resulting from {it:c(level)} (see {helpb level}).

{phang} 
{opt rseed(#)} sets the random-number seed to ensure reproducibility of the prediction intervals. 

{marker examples}{...}
{title:Examples}

{hline}
{pstd}
In-sample statistics{p_end}
{hline}

{pstd}
Load data and estimation
{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. hetwals price (c.mpg) c.mpg#c.mpg length i.foreign, het(length i.foreign c.mpg##c.mpg)}{p_end}

{pstd}
Linear prediction
{p_end}
{phang2}{cmd:. predict double hetwals_lp}{p_end}
{phang2}{cmd:. sum hetwals_lp}{p_end}

{pstd}
Estimated standard deviations of the errors{p_end}
{phang2}{cmd:. predict double hetwals_sd, sigma}{p_end}
{phang2}{cmd:. sum hetwals_sd}{p_end}

{pstd}
90% Confidence intervals for E(y_i,x_i){p_end}
{phang2}{cmd:. predict double hetwals_ci_lb hetwals_ci_ub, cint level(90)}{p_end}
{phang2}{cmd:. sum hetwals_ci_lb hetwals_ci_ub}{p_end}

{hline}
{pstd}
Out-of-sample statistics {p_end}
{hline}

{pstd}
Load data and estimation{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. hetwals price (c.mpg) c.mpg#c.mpg length i.foreign if _n!=1, het(length i.foreign c.mpg##c.mpg)}{p_end}

{pstd}
Linear prediction{p_end}
{phang2}{cmd:. predict double hetwals_lp_2 if e(sample)!=1}{p_end}
{phang2}{cmd:. sum hetwals_lp_2}{p_end}

{pstd}
Estimated standard deviations of the errors{p_end}
{phang2}{cmd:. predict double hetwals_sd_2 if e(sample)!=1, sigma}{p_end}
{phang2}{cmd:. sum hetwals_sd_2}{p_end}

{pstd}
90% Prediction intervals for y_i|(x_i,v_i){p_end}
{phang2}{cmd:. predict double hetwals_pi_lb hetwals_pi_ub if e(sample)!=1, pint rseed(12345) level(90)}{p_end}
{phang2}{cmd:. list price hetwals_lp_2 hetwals_sd_2 hetwals_pi_lb hetwals_pi_ub if _n==1}{p_end}
 
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

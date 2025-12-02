{smcl}
{* *! version 1.0.0  04oct2024}{...}
{vieweralsosee "glmwals" "help glmwals"}{...}
{vieweralsosee "glmwals postestimation" "help glmwals_postestimation"}{...}
{vieweralsosee "lcwals" "help lcwals"}{...}
{vieweralsosee "margwals" "help margwals"}{...}
{vieweralsosee "wals" "help wals"}{...}
{vieweralsosee "hetwals" "help hetwals"}{...}
{vieweralsosee "ar1wals" "help ar1wals"}{...}
{vieweralsosee "xtwals" "help xtwals"}{...}
{viewerjumpto "Syntax" "glmwals_p##syntax"}{...}
{viewerjumpto "Description" "glmwals_p##description"}{...}
{viewerjumpto "Examples" "glmwals_p##examples"}{...}
{viewerjumpto "References" "glmwals_p##references"}{...}
{p2colset 1 12 14 2}{...}
{p2col:{bf:predict} {hline 2}}Weighted-average least squares (WALS) predictions for generalized linear models (GLMs){p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:predict} {dtype} {newvar} {ifin} 
[{cmd:,} {help glmwals_p##opts1:{it:predict_options}}]

{synoptset 20 tabbed}{...}
{marker opts1}{...}
{synopthdr :predict_options}
{synoptline}
{synopt :{opt mu}}plug-in estimate of the expected value of {depvar} (default){p_end}
{synopt :{opt xb}}linear prediction{p_end}
{synopt :{opt biasp}}bias of the linear prediction{p_end}
{synopt :{opt stdp}}standard error of the linear prediction{p_end}
{synopt :{opt rmsep}}root mean squared error of the linear prediction{p_end}
{synopt :{opt bcp}}bias-corrected linear prediction{p_end}
{synopt :{opt stdbcp}}standard error of the bias-corrected linear prediction{p_end}
{synopt :{opt r:esiduals}}difference between observed and fitted outcomes{p_end}
{synopt :{opt nooff:set}}modify calculations to ignore offset variable{p_end}
{synoptline}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd}
{cmd:predict} uses the current estimation results produced by {helpb glmwals} to compute weighted-average least squares (WALS) predictions 
such as expected values, linear predictions, bias-corrected linear predictions, sampling moments (bias, SE, and RMSE) of linear pridictions
and the bias-corrected linear predictions, 
and residuals. 

{pstd}
By default all statistics are calculated out of sample; type {cmd:predict ... if e(sample) ...} 
if calculations must be restricted to the estimation sample.
For additional details see {help glmwals_p##delucamagnus2024a:De Luca and Magnus (2024a,} {help glmwals_p##delucamagnus2024b:2024b)}.
To compute point estimates and simulation-based confidence intervals for mean responses and marginal effects
at given values of the regressors see also {help margwals}.

{marker options}{...}
{title:Options}

{dlgtab:predict_options1}

{phang} 
{opt mu}, the default, calculates the plug-in estimate of the expected value of {depvar}.

{phang} 
{opt xb} calculates the linear prediction.

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
{opt r:esiduals} calculates the so-called response residuals, that is the difference between the observed outcome and its expected value.

{phang} 
{opt nooff:set} modifies the default calculations made by {com: predict} to ignore the offset variable (if applicable).

{marker examples}{...}
{title:Examples}

{hline}
{pstd}
Binomial model with logit link {p_end}
{hline}

{pstd}
Setup and estimation
{p_end}
{phang2}{cmd:. webuse lbw, clear}{p_end}
{phang2}{cmd:. glmwals low (lwt) age i.race smoke ptl ht ui, family(binomial) prior(cauchy) rseed(12345)}{p_end}

{pstd}
Plug-in estimate of the response probability
{p_end}
{phang2}{cmd:. predict glmwals_mu}{p_end}

{pstd}
Linear prediction
{p_end}
{phang2}{cmd:. predict glmwals_xb, xb}{p_end}

{pstd}
Plug-in estimate of the bias of the linear predictor
{p_end}
{phang2}{cmd:. predict glmwals_biasp, biasp}{p_end}

{pstd}
Plug-in estimate of the SE of the linear predictor
{p_end}
{phang2}{cmd:. predict glmwals_stdp, stdp}{p_end}

{pstd}
Plug-in estimate of the RMSE of the linear predictor
{p_end}
{phang2}{cmd:. predict glmwals_rmsep, rmsep}{p_end}



{hline}
{pstd}
Poisson model with log link {p_end}
{hline}

{pstd}
Setup and estimation
{p_end}
{phang2}{cmd:. use "https://www.stata-press.com/data/r18/breathe", clear}{p_end}
{phang2}{cmd:. do https://www.stata-press.com/data/r18/no2}{p_end}
{phang2}{cmd:. display "$cc"}{p_end}
{phang2}{cmd:. display "$fc"}{p_end}
{phang2}{cmd:. rename no2_class NO2}{p_end}
{phang2}{cmd:. glmwals omissions (NO2) $cc i.($fc), family(poisson) eform noaux rseed(12345)}{p_end}
 
{pstd}
Predicted mean
{p_end}
{phang2}{cmd:. predict double mean_om}{p_end}


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

{smcl}
{* *! version 1.0.0  04oct2024}{...}
{vieweralsosee "wals" "help wals"}{...}
{vieweralsosee "wals postestimation" "help wals_postestimation"}{...}
{vieweralsosee "lcwals" "help lcwals"}{...}
{vieweralsosee "margwals" "help margwals"}{...}
{vieweralsosee "hetwals" "help hetwals"}{...}
{vieweralsosee "ar1wals" "help ar1wals"}{...}
{vieweralsosee "xtwals" "help xtwals"}{...}
{vieweralsosee "glmwals" "help glmwals"}{...}
{viewerjumpto "Syntax" "wals_p##syntax"}{...}
{viewerjumpto "Description" "wals_p##description"}{...}
{viewerjumpto "Examples" "wals_p##examples"}{...}
{viewerjumpto "References" "wals_p##references"}{...}
{p2colset 1 12 14 2}{...}
{p2col:{bf:predict} {hline 2}}Weighted-average least squares (WALS) predictions for linear models with i.i.d. errors{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Syntax of the {com:predict} command for {helpb wals} is presented under the following headings:

{phang2}{help wals_p##syntax1:Compute linear predictions, moments of linear predictions, and residuals}{p_end}
{phang2}{help wals_p##syntax2:Compute confidence and prediction intervals}{p_end}


{marker syntax1}{...}
{phang}{ul:{bf:Compute linear predictions, moments of linear predictions, and residuals}}

{p 8 16 2}
{cmd:predict} {dtype} {newvar} {ifin} 
[{cmd:,} {help wals_p##opts1:{it:predict_options1}}]

{marker syntax2}{...}
{phang}{ul:{bf:Compute confidence and prediction intervals}}

{p 8 16 2}
{cmd:predict} {dtype} {newvar}_{it:l}  {newvar}_{it:u} {ifin} 
[{cmd:,} {help wals_p##opts2:{it:predict_options2}}]


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
{cmd:predict} uses the current estimation results produced by {helpb wals} to compute weighted-average least squares (WALS) predictions 
such as linear predictions, bias-corrected predictions, their sampling moments (bias, SE, and RMSE), residuals, bias-corrected residuals,
confidence intervals for E(y_i|x_i), and prediction intervals for y_i|x_i. 

{pstd}
In particular, the first syntax creates in {newvar} one of the statistics available in {help wals_p##opts1:{it:predict_options1}}.
If no option is specified, the default {com: xb} option is assumed. 
Estimated sampling moments and bias-corrected statistics depend on the plug-in estimators specified in the last {helpb wals} regression 
(see {help wals##options:{it:wals options}}).

{pstd}
The second syntax creates in {newvar}_{it:l} and {newvar}_{it:u} 
the lower and upper bounds of the confidence intervals for E(y_i|x_i) if one specifies the {cmd:cinterval} option,
and the lower and upper bounds of the prediction intervals for y_i|x_i if one specifies the {cmd:pinterval} option.
The confidence level of both intervals is controlled by the {opt level(#)} option.

{pstd}
Except for the linear prediction, all statistics require that 
the last {helpb wals} regression was fitted without the {opt fast} option (see {help wals##options:{it:wals options}}).
All statistics are available both in sample and out of sample (default); type {cmd:predict ... if e(sample) ...} 
if calculations must be restricted to the estimation sample.
For additional details see {help wals_p##delucamagnus2024:De Luca and Magnus (2024)}.

{marker options}{...}
{title:Options}

{dlgtab:predict_options1}

{phang} 
{opt xb}, the default, calculates the linear prediction (i.e., the WALS point estimate of E(y_i|x_i)).

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

{dlgtab:predict_options2}

{phang} 
{opt cint:erval} calculates the confidence intervals for E(y_i|x_i).

{phang} 
{opt pint:erval} calculates the prediction intervals for y_i|x_i. 

{phang} 
{opt lev:el(#)} specifies the confidence level of the confidence/prediction intervals, where # is a percentage between 10 and 99.99.
The default is {com:level(95)} or the confidence level resulting from {it:c(level)} (see {helpb level}).

{phang} 
{opt rseed(#)} sets the random-number seed to ensure reproducibility of the prediction intervals. 

{marker examples}{...}
{title:Examples}

{hline}
{pstd}
In-sample predictions based on plug-in double-shrinkage estimates of the sampling moments{p_end}
{hline}

{pstd}
Setup and estimation{p_end}
{phang2}{cmd:. use "https://www.stata-press.com/data/r18/uscrime", clear}{p_end}
{phang2}{cmd:. wals ln_offenses ln_malepop-ln_prisont, rseed(1234)}{p_end}

{pstd}
Linear prediction{p_end}
{phang2}{cmd:. predict double wals_lp_1}{p_end}
{phang2}{cmd:. sum wals_lp_1}{p_end}

{pstd}
Bias-corrected linear prediction{p_end}
{phang2}{cmd:. predict double wals_bclp_1, bcp}{p_end}
{phang2}{cmd:. sum wals_bclp_1}{p_end}

{pstd}
Residuals from the linear prediction{p_end}
{phang2}{cmd:. predict double wals_res_1, res}{p_end}
{phang2}{cmd:. sum wals_res_1}{p_end}

{pstd}
Residuals from the bias-corrected linear prediction{p_end}
{phang2}{cmd:. predict double wals_bcr_1, bcres}{p_end}
{phang2}{cmd:. sum wals_bcr_1}{p_end}

{pstd}
99% confidence intervals for E(y_i|x_i){p_end}
{phang2}{cmd:. predict double wals_lci_1 wals_uci_1, cint level(99)}{p_end}
{phang2}{cmd:. format %9.3f ln_offenses wals_*_1}{p_end}
{phang2}{cmd:. ln_offenses wals_lci_1 wals_lp_1 wals_bclp_1 wals_uci_1 wals_res_1 wals_bcr_1 in 1/5,abbr(12)}{p_end}

{hline}
{pstd}
Out-of-sample predictions based on plug-in maximum likelihood estimates of the sampling moments{p_end}
{hline}

{pstd}
Setup and estimation{p_end}
{phang2}{cmd:. splitsample, generate(sample) nsplit(2) split(2 1) rseed(18)}{p_end}
{phang2}{cmd:. wals ln_offenses (ln_malepop) southern-ln_prisont if sample==1, plugin(ml) rseed(1234)}{p_end}

{pstd}
Linear prediction{p_end}
{phang2}{cmd:. predict double wals_lp_2 if sample==2}{p_end}
{phang2}{cmd:. sum wals_lp_2}{p_end}

{pstd}
Bias of the linear prediction{p_end}
{phang2}{cmd:. predict double wals_bias_lp_2 if sample==2, biasp}{p_end}
{phang2}{cmd:. sum wals_bias_lp_2}{p_end}

{pstd}
SE of the linear prediction{p_end}
{phang2}{cmd:. predict double wals_se_lp_2 if sample==2, stdp}{p_end}
{phang2}{cmd:. sum wals_se_lp_2}{p_end}

{pstd}
RMSE of the linear prediction{p_end}
{phang2}{cmd:. predict double wals_rmse_lp_2 if sample==2, rmsep}{p_end}
{phang2}{cmd:. sum wals_rmse_lp_2}{p_end}

{pstd}
SE of the bias-corrected linear prediction{p_end}
{phang2}{cmd:. predict double wals_se_bclp_2 if sample==2, stdbcp}{p_end}
{phang2}{cmd:. sum wals_se_bclp_2}{p_end}

{pstd}
99% prediction interval for y_i|x_i{p_end}
{phang2}{cmd:. predict double wals_lpi_2 wals_upi_2 if sample==2, pint level(99) rseed(1234)}{p_end}
{phang2}{cmd:. format %9.3f ln_offenses wals_*_2}{p_end}
{phang2}{cmd:. list ln_offenses wals_lpi_2 wals_lp_2 wals_upi_2 if sample==2,abbr(12)}{p_end} 
 
{marker references}{...}
{title:References}

{marker delucamagnus2024}{...}
{phang}
De Luca, G., and J.R. Magnus. 2024.
Weighted-Average Least Squares: Improvements and extensions. 
{it:Stata Journal}, forthcoming.
Available at: 
{browse "https://www.janmagnus.nl/"}

{title:Authors}

{pstd}Giuseppe De Luca{p_end}
{pstd}University of Palermo{p_end}
{pstd}Palermo, Italy{p_end}
{pstd}emal:giuseppe.deluca@unipa.it{p_end}

{pstd}Jan R. Magnus {p_end}
{pstd}Vrije Universiteit Amsterdam{p_end}
{pstd}Amsterdam, The Netherlands{p_end}
{pstd}emal:jan@janmagnus.nl{p_end}


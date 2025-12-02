{smcl}
{* *! version 1.0.0  04oct2024}{...}
{vieweralsosee "ar1wals" "help ar1wals"}{...}
{vieweralsosee "ar1wals postestimation" "help ar1wals_postestimation"}{...}
{vieweralsosee "lcwals" "help lcwals"}{...}
{vieweralsosee "margwals" "help margwals"}{...}
{vieweralsosee "wals" "help wals"}{...}
{vieweralsosee "hetwals" "help hetwals"}{...}
{vieweralsosee "xtwals" "help xtwals"}{...}
{vieweralsosee "glmwals" "help glmwals"}{...}
{viewerjumpto "Syntax" "ar1wals_p##syntax"}{...}
{viewerjumpto "Description" "ar1wals_p##description"}{...}
{viewerjumpto "Examples" "ar1wals_p##examples"}{...}
{viewerjumpto "References" "ar1wals_p##references"}{...}
{p2colset 1 12 14 2}{...}
{p2col:{bf:predict} {hline 2}}Weighted-average least squares (WALS) predictions for linear models with stationary AR(1) errors{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Syntax of the {com:predict} command for {helpb ar1wals} is presented under the following headings:

{phang2}{help ar1wals_p##syntax1:Compute linear predictions, moments of linear predictions, and residuals}{p_end}
{phang2}{help ar1wals_p##syntax2:Compute confidence and prediction intervals}{p_end}


{marker syntax1}{...}
{phang}{ul:{bf:Compute linear predictions, moments of linear predictions, and residuals}}

{p 8 16 2}
{cmd:predict} {dtype} {newvar} {ifin} 
[{cmd:,} {help ar1wals_p##opts1:{it:predict_options1}}]

{marker syntax2}{...}
{phang}{ul:{bf:Compute confidence and prediction intervals}}

{p 8 16 2}
{cmd:predict} {dtype} {newvar}_{it:l}  {newvar}_{it:u} {ifin} 
[{cmd:,} {help ar1wals_p##opts2:{it:predict_options2}}]


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
{synopt :{opt wp}}{it:s} periods-ahead WALS predictions{p_end}
{synopt :{opt bcwp}}{it:s} periods-ahead bias-corrected WALS predictions{p_end}
{synoptline}
{p2colreset}{...}


{synoptset 20 tabbed}{...}
{marker opts2}{...}
{synopthdr :predict_options2}
{synoptline}
{synopt :{opt cint:erval}}confidence intervals for E(y_{T+s}|x_{T+s}) {p_end}
{synopt :{opt pint:erval}}prediction intervals for y_{T+s}|x_{T+s}{p_end}
{synopt :{opt lev:el(#)}}sets the confidence level{p_end}
{synopt :{opt rseed(#)}}sets the random-number seed for the prediction intervals{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:predict} uses the current estimation results produced by {helpb ar1wals} to compute weighted-average least squares (WALS) predictions 
such as linear predictions, bias-corrected predictions, their sampling moments (bias, SE, and RMSE), residuals, bias-corrected residuals,
{it:s} periods-ahead WALS predictions, confidence intervals for E(y_{T+s}|x_{T+s}), and prediction intervals for y_{T+s}|x_{T+s}. 

{pstd}
In particular, the first syntax creates in {newvar} one of the statistics available in {help ar1wals_p##opts1:{it:predict_options1}}.
If no option is specified, the default {com: xb} option is assumed. 
Estimated sampling moments and bias-corrected statistics depend on the plug-in estimators specified in the last {helpb ar1wals} regression 
(see {help wals##options:{it:wals_options}}).
All statistics are conditional on the first-step estimate of the autocorrelation coefficient ρ and
the only statistic allowed with the {help margwals} command is the default linear prediction.

{pstd}
The second syntax creates in {newvar}_{it:l} and {newvar}_{it:u} 
the lower and upper bounds of the confidence intervals for E(y_{T+s}|x_{T+s}) if one specifies the {cmd:cinterval} option,
and the lower and upper bounds of the prediction intervals for y_{T+s}|x_{T+s} if one specifies the {cmd:pinterval} option.
The confidence level of both intervals is controlled by the {opt level(#)} option.

{pstd}
Except for the {opt xb} and {opt wp} options, all statistics require that 
the last {helpb ar1wals} regression was fitted without the {opt fast} option (see {help wals##options:{it:wals_options}}).
All statistics are available both in sample and out of sample (default); type {cmd:predict ... if e(sample) ...} 
if calculations must be restricted to the estimation sample.
For additional details see {help ar1wals_p##delucamagnus2024a:De Luca and Magnus (2024a,} {help ar1wals_p##delucamagnus2024b:2024b)}.

{marker options}{...}
{title:Options}

{dlgtab:predict_options1}

{phang} 
{opt xb}, the default, calculates the linear prediction without accounting for the estimated correlation among residuals.

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
{opt wp} calculates the {it:s} periods-ahead WALS predictions accounting for the estimated correlation among residuals 
as in {help ar1wals_p##goldberger1962:Goldberger (1962)}.

{phang} 
{opt bcwp} calculates the {it:s} periods-ahead bias-corrected WALS predictions accounting for the estimated correlation among residuals 
as in {help ar1wals_p##goldberger1962:Goldberger (1962)}.

{dlgtab:predict_options2}

{phang} 
{opt cint:erval} calculates the confidence intervals for E(y_{T+s}|x_{T+s}).

{phang} 
{opt pint:erval} calculates the prediction intervals for y_{T+s}|x_{T+s}. 

{phang} 
{opt lev:el(#)} specifies the confidence level of the confidence/prediction intervals, where # is a percentage between 10 and 99.99.
The default is {com:level(95)} or the confidence level resulting from {it:c(level)} (see {helpb level}).

{phang} 
{opt rseed(#)} sets the random-number seed to ensure reproducibility of the prediction intervals. 

{marker examples}{...}
{title:Examples}

{hline}
{pstd}
barium data from {help ar1wals_p##wooldridge2013:Wooldridge (2013, Example 10.11)} {p_end}
{hline}

{pstd}
Setup
{p_end}
{phang2}{cmd:. bcuse barium, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. sort t}{p_end}
{phang2}{cmd:. gen esample=(_n<=_N-24)}{p_end}

{pstd}
LS predictions
{p_end}
{phang2}{cmd:. regress lchnimp lchempi lrtwex lgas befile6 affile6 afdec6 feb-dec if esample==1}{p_end}
{phang2}{cmd:. predict double ols_p if e(sample)!=1}{p_end}

{pstd}
Best linear unbiased predictions ({help ar1wals_p##goldberger1962:Goldberger 1962})
{p_end}
{phang2}{cmd:. prais lchnimp lchempi lrtwex lgas befile6 affile6 afdec6 feb-dec if esample==1}{p_end}
{phang2}{cmd:. predict double prais_xb if e(sample)!=1}{p_end}
{phang2}{cmd:. sort t}{p_end}
{phang2}{cmd:. predict double prais_res if esample==1 & _n==e(N), res}{p_end}
{phang2}{cmd:. sum t if esample==1 & _n==e(N), meanonly}{p_end}
{phang2}{cmd:. gen pow_s=t-r(mean) if e(sample)!=1}{p_end}
{phang2}{cmd:. sum prais_res, meanonly}{p_end}
{phang2}{cmd:. gen double prais_p=prais_xb+e(rho)^pow_s * r(mean) if e(sample)!=1}{p_end}

{pstd}
WALS predictions
{p_end}
{phang2}{cmd:. ar1wals lchnimp (lchempi lrtwex lgas) befile6 affile6 afdec6 feb-dec if esample==1, rseed(12345) prior(log)}{p_end}
{phang2}{cmd:. predict double ar1wals_p if e(sample)!=1, wp}{p_end}

{pstd}
Mean squared prediction errors
{p_end}
{phang2}{cmd:. gen double mspe_ols=(ols_p-lchnimp)^2 if esample!=1}{p_end}
{phang2}{cmd:. gen double mspe_prais=(prais_p-lchnimp)^2 if esample!=1}{p_end}
{phang2}{cmd:. gen double mspe_ar1wals=(ar1wals_p-lchnimp)^2 if esample!=1}{p_end}
{phang2}{cmd:. tabstat mspe_*, stat(mean) col(stat) varwidth(15)}{p_end}

{pstd}
95% prediction intervals
{p_end}
{phang2}{cmd:. predict double ar1wals_lb ar1wals_ub , pinterval rseed(12345)}{p_end}
 
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

{marker goldberger1962}{...}
{phang}
Goldberger, A. S. 1962.
Best linear unbiased prediction in the generalized linear regression model.
{it:Journal of the American Statistical Association} 57: 369–375.
{browse "https://doi.org/10.2307/2281645"}

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

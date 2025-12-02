{smcl}
{* *! version 1.0.0  04oct2024}{...}
{vieweralsosee "wals" "help wals"}{...}
{vieweralsosee "wals postestimation" "help wals_postestimation"}{...}
{vieweralsosee "predict" "help wals_p"}{...}
{vieweralsosee "margwals" "help margwals"}{...}
{vieweralsosee "hetwals" "help hetwals"}{...}
{vieweralsosee "ar1wals" "help ar1wals"}{...}
{vieweralsosee "xtwals" "help xtwals"}{...}
{vieweralsosee "glmwals" "help glmwals"}{...}
{viewerjumpto "Syntax" "lcwals##syntax"}{...}
{viewerjumpto "Description" "lcwals##description"}{...}
{viewerjumpto "Options" "lcwals##options"}{...}
{viewerjumpto "Examples" "lcwals##examples"}{...}
{viewerjumpto "Stored results" "lcwals##results"}{...}
{viewerjumpto "References" "lcwals##references"}{...}
{vieweralsosee "[R] lincom" "help lincom"}{...}
{p2colset 1 11 13 2}{...}
{p2col:{bf:lcwals} {hline 2}}Linear combination of weighted-average least squares (WALS) estimates{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:lcwals} {it:{help exp}} [{cmd:,} {it:options}]

{synoptset 16}{...}
{synopthdr}
{synoptline}
{synopt :{opt level(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt :{opt ef:orm}}reports the results in exponential form{p_end}
{synopt :{opt lones:ide}}the left one-sided confidence intervals; default is two-sided confidence intervals{p_end}
{synopt :{opt rones:ide}}the right one-sided confidence intervals; default is two-sided confidence intervals{p_end}
{synopt :{opt cformat(%fmt)}}set the display format for statistics in the coefficient table; default is %8.0g{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}{it:exp} is any linear combination of coefficients that is a valid
syntax for {cmd:lincom}; see {helpb lincom:[R] lincom}. 
However, arithmetical operators are restricted to +, -, *, and /,
parentheses are not allowed, 
and coefficients of the linear combination must be specified as numerical values without using
other mathematical functions.
{it:exp} must not contain an equal sign.

{marker description}{...}
{title:Description}

{pstd}
{cmd:lcwals} uses the latest estimation results from {helpb wals} to compute the weighted-average least squares (WALS) estimate
of a linear combination of the focus and auxiliary parameters.
In addition to the point estimate, it calculates estimates of the sampling moments (bias, SE, and RMSE) and simulation-based confidence 
intervals for the required linear combinations of the focus and auxiliary parameters. 

{pstd}
Specifying the {com:loneside} and {com: roneside} options, one can also compute the left and right one-sided
confidence intervals instead of the (default) two-sided confidence intervals.
See {help lcwals##delucamagnus2024a:De Luca and Magnus (2024a)} for additional details.

{pstd}
{cmd:lcwals} can also be used after the {helpb hetwals}, {helpb ar1wals}, {helpb xtwals}, and {helpb glmwals} commands
discussed in {help lcwals##delucamagnus2024b:De Luca and Magnus (2024b)}. 
 
{marker options}{...}
{title:Options}

{phang} 
{opt level(#)} specifies the confidence level, as a percentage, for confidence intervals.
The default is {com:level(95)} or as set by {com:set level}.

{phang} 
{opt eform} reports coefficient estimates in exponential form, that is as exp(b) rather than b.
Estimated sampling moments and confidence intervals are similarly transformed.
This option is seldom used in the context of linear regression models, but it could be useful 
in other classes of models such as GLM. 

{phang} 
{opt loneside} compute the left one-sided confidence intervals instead of the default two-sided confidence
intervals.

{phang} 
{opt roneside} compute the right one-sided confidence intervals instead of the default two-sided confidence
intervals.

{phang}
{opt cformat(%fmt)} specifies the display format for coefficients, estimated moments, and confidence limits 
in the coefficient table. 
The default format is %8.0g and the maximum format width is 8. 

{marker examples}{...}
{title:Examples}

{hline}
{pstd}
OLS estimates {p_end}
{hline}

{pstd}
Setup and estimation{p_end}
{phang2}{cmd:. webuse regress, clear}{p_end}
{phang2}{cmd:. regress y x1 x2 x3}{p_end}

{pstd}
OLS estimates of linear combiantions of the regression parameters{p_end}
{phang2}{cmd:. lincom x2-x1}{p_end}
{phang2}{cmd:. lincom 3*x1 + 500*x3}{p_end}
{phang2}{cmd:. lincom 3*x1 + 500*x3 - 12, level(99)}{p_end}

{hline}
{pstd}
WALS estimates {p_end}
{hline}

{pstd}
WALS estimates of the regression parameters{p_end}
{phang2}{cmd:. wals y x1 x2 (x3), rseed(1234)}{p_end}

{pstd}
WALS estimates of linear combiantions of the regression parameters{p_end}
{phang2}{cmd:. lcwals x2-x1}{p_end}
{phang2}{cmd:. lcwals 3*x1 + 500*x3}{p_end}
{phang2}{cmd:. lcwals 3*x1 + 500*x3 - 12, level(99)}{p_end}
{phang2}{cmd:. return list}{p_end}
{phang2}{cmd:. lcwals 3*x1 + 500*x3 - 12, loneside level(99)}{p_end}
{phang2}{cmd:. return list}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:lcwals} stores the following in {cmd:r()}:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}

{synopt:{cmd:r(estimate)}}point estimate of the linear combination{p_end}
{synopt:{cmd:r(bias)}}estimated bias{p_end}
{synopt:{cmd:r(se)}}estimated SE{p_end}
{synopt:{cmd:r(rmse)}}estimated RMSE{p_end}
{synopt:{cmd:r(lb)}}lower bound of confidence interval; . if {com:loneside} is specified{p_end}
{synopt:{cmd:r(ub)}}upper bound of confidence interval; . if {com:roneside} is specified{p_end}
{synopt:{cmd:r(level)}}confidence level{p_end}
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

{title:Authors}

{pstd}Giuseppe De Luca{p_end}
{pstd}University of Palermo{p_end}
{pstd}Palermo, Italy{p_end}
{pstd}emal:giuseppe.deluca@unipa.it{p_end}

{pstd}Jan R. Magnus {p_end}
{pstd}Vrije Universiteit Amsterdam{p_end}
{pstd}Amsterdam, The Netherlands{p_end}
{pstd}emal:jan@janmagnus.nl{p_end}

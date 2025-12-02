{smcl}
{* *! version 1.0.0  04oct2024}{...}
{vieweralsosee "ar1wals" "help ar1wals"}{...}
{vieweralsosee "predict" "help ar1wals_p"}{...}
{vieweralsosee "lcwals" "help lcwals"}{...}
{vieweralsosee "margwals" "help margwals"}{...}
{vieweralsosee "wals" "help wals"}{...}
{vieweralsosee "hetwals" "help hetwals"}{...}
{vieweralsosee "xtwals" "help xtwals"}{...}
{vieweralsosee "glmwals" "help glmwals"}{...}
{vieweralsosee "[R] estimates" "help estimates"}{...}
{p2colset 1 28 30 2}{...}
{p2col:{bf:ar1wals postestimation} {hline 2}}Postestimation tools for 
weighted-average least squares (WALS) estimation of linear models with stationary AR(1) errors{p_end}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd}
The following postestimation commands are available after {helpb ar1wals}:

{synoptset 24 tabbed}{...}
{p2coldent :Command}Description{p_end}
{synoptline}
{synopt :{helpb ar1wals_p:predict}}WALS predictions{p_end}
{synopt :{helpb lcwals:lcwals}}Linear combinations of WALS estimates{p_end}
{synopt :{helpb margwals:margwals}}WALS estimates of marginal means, predictive margins, marginal effects, and average marginal effects{p_end}
{p2coldent:* {helpb estimates}}cataloging estimation results{p_end}
{synoptline}
{p 4 6 2}
* {cmd:estimates stats} is not appropriate after {cmd:ar1wals} estimation.

{pstd}
The above sections are not included in this help file.
{p_end}

{title:Authors}

{pstd}Giuseppe De Luca{p_end}
{pstd}University of Palermo{p_end}
{pstd}Palermo, Italy{p_end}
{pstd}emal:giuseppe.deluca@unipa.it{p_end}

{pstd}Jan R. Magnus {p_end}
{pstd}Vrije Universiteit Amsterdam{p_end}
{pstd}Amsterdam, The Netherlands{p_end}
{pstd}emal:jan@janmagnus.nl{p_end}

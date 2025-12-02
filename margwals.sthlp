{smcl}
{* *! version 1.0.0  04oct2024}{...}
{vieweralsosee "wals" "help wals"}{...}
{vieweralsosee "wals postestimation" "help wals_postestimation"}{...}
{vieweralsosee "predict" "help wals_p"}{...}
{vieweralsosee "lcwals" "help lcwals"}{...}
{vieweralsosee "hetwals" "help hetwals"}{...}
{vieweralsosee "ar1wals" "help ar1wals"}{...}
{vieweralsosee "xtwals" "help xtwals"}{...}
{vieweralsosee "glmwals" "help glmwals"}{...}
{vieweralsosee "[R] margins" "help margins"}{...}
{viewerjumpto "Syntax" "margwals##syntax"}{...}
{viewerjumpto "Description" "margwals##description"}{...}
{viewerjumpto "Options" "margwals##options"}{...}
{viewerjumpto "Examples" "margwals##examples"}{...}
{viewerjumpto "Stored results" "margwals##results"}{...}
{viewerjumpto "References" "lcwals##references"}{...}
{p2colset 1 13 15 2}{...}
{p2col:{bf:margwals} {hline 2}}Weighted-average least squares (WALS) estimates of marginal means, predictive margins, marginal effects, and average marginal effects{p_end}
{p2colreset}{...}


{marker syntax_margwals}{...}
{marker margwals}{...}
{title:Syntax}

{p 8 16 2}
{cmd:margwals} [{help fvvarlist:{it:marginlist}}] 
{ifin} 
[{it:{help margwals##weight:weight}}]
[{cmd:, }{help margwals##optstbl:{it:margwals_options}}]

{p 4 6 2}
where {it:marginlist} is a list of factor variables or interactions that appear in the latest
WALS regression, and {it:margwals_options} are the following options

{marker margwals_options}{...}
{synoptset 22 tabbed}{...}
{synopthdr:margwals_options}
{synoptline}
{syntab :Response options}
{synopt:{opt pr:edict(pred_opt)}}estimate
	margins for {cmd:predict,} {it:pred_opt}{p_end}
{synopt:{opth dydx(varlist)}}estimate
	marginal effect of variables in {it:varlist}{p_end}
{synopt:{opth eyex(varlist)}}estimate
	elasticities of variables in {it:varlist}{p_end}
{synopt:{opth dyex(varlist)}}estimate
	semielasticity -- d(y)/d(ln x){p_end}
{synopt:{opth eydx(varlist)}}estimate
	semielasticity -- d(ln y)/d(x){p_end}
{synopt:{opt cont:inuous}}treat factor-level indicators as continuous{p_end}

{syntab :Other options from the {helpb margins} command}
{synopt:{opt grand}}add 
	the overall margin; default if no {it:marginlist}{p_end}
{synopt:{cmd:at(}{it:{help margins##atspec:atspec}{cmd:)}}}estimate 
	margins at specified values of covariates{p_end}
{synopt:{opt atmeans}}estimate 
	margins at the means of covariates{p_end}
{synopt:{opt asbal:anced}}treat
	all factor variables as balanced{p_end}
{synopt:{opth over(varlist)}}estimate
	margins at unique values of {it:varlist}{p_end}
{synopt:{opth within(varlist)}}estimate
	margins at unique values of the nesting factors in {it:varlist}{p_end}
{synopt:{opt noweight:s}}ignore
	weights specified in estimation{p_end}
{synopt:{opt noe:sample}}do
	not restrict {cmd:margins} to the estimation sample{p_end}
{synopt:{opt emptycells}{cmd:(}{it:{help margins##empspec:empspec}{cmd:)}}}treatment of empty cells for balanced factors{p_end}
{synopt:{opt estimtol:erance(tol)}}specify numerical tolerance used to determine estimable functions; default is {cmd:estimtolerance(1e-5)}{p_end}
{synopt:{opt noestimcheck}}suppress estimability
	checks{p_end}
{synopt :{opt force}}estimate margins
	despite potential problems{p_end}
{synopt :{opt chain:rule}}use the
	chain rule when computing derivatives{p_end}
{synopt :{opt nochain:rule}}do
	not use the chain rule{p_end}
{synopt:{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:margwals} uses the latest estimation results from {helpb wals} to compute the weighted-average least squares (WALS) estimates of 
smooth, possibly nonlinear, real-valued functions g(β,x) of the parameters β at some value x of the regressors.
In addition to point estimates, it provides a simulation-based confidence interval for g(β,x) 
using the Monte Carlo replications of the latest bias-corrected WALS estimates of β.
See {help margwals##delucamagnus2024a:De Luca and Magnus (2024a)} for additional details.
  
{pstd}
As for the {helpb margins} command, capabilities include marginal means, predictive margins,
marginal effects, and average marginal effects.
Estimates of sampling moments, in particular the bias, of the WALS estimator of g(β,x) are not currently available. 
For WALS estimation of linear combinations of the focus and auxiliary parameters see {helpb lcwals}.
 
{pstd}
{cmd:margwals} can also be used after the {helpb hetwals}, {helpb ar1wals}, {helpb xtwals}, and {helpb glmwals} commands
discussed in {help margwals##delucamagnus2024b:De Luca and Magnus (2024b)}. 

{marker options}{...}
{title:Options}

{pstd}
For a detailed description of the various options available in {com:margwals} see the corresponding {help margins##options:options of margins}.
Note that some options of {com:margins} have been excluded from {com:margwals} either because they may lead to a misleading interpretation of the estimation results, or because they cannot be 
currently applied to the WALS estimator.  

{marker examples}{...}
{title:Examples}

{pstd}WALS margins of responses{p_end}
{phang2}{cmd:. webuse margex, clear}{p_end}
{phang2}{cmd:. regress y i.sex i.group i.sex#i.group}{p_end}
{phang2}{cmd:. margins sex, level(99)}{p_end}
{phang2}{cmd:. wals y (i.sex i.group) i.sex#i.group, rseed(12345)}{p_end}
{phang2}{cmd:. margins sex, level(99)        // invalid p-values and confidence intervals}{p_end}
{phang2}{cmd:. margwals sex, level(99)      // simulation-based confidence intervals}{p_end}

{pstd}Multiple WALS margins from one {com:margwals} command{p_end}
{phang2}{cmd:. margwals sex group}{p_end}

{pstd}WALS margins with continuous variables{p_end}
{phang2}{cmd:. wals y (i.sex i.group) i.sex#i.group age, rseed(12345)}{p_end}
{phang2}{cmd:. margwals sex group}{p_end}
{phang2}{cmd:. margwals, at(age=(30 35 40 45 50))}{p_end}

{pstd}WALS margins of interactions{p_end}
{phang2}{cmd:. margwals sex#group}{p_end}

{pstd}WALS average marginal effects of all covariates{p_end}
{phang2}{cmd:. margwals, dydx(*)}{p_end}

{pstd}Obtaining WALS margins as though the data were balanced{p_end}
{phang2}{cmd:. webuse acmemanuf, clear}{p_end}
{phang2}{cmd:. regress y pressure##temp}{p_end}
{phang2}{cmd:. margins, asbalanced}{p_end}
{phang2}{cmd:. wals y (i.pressure i.temp) pressure#temp, rseed(12345)}{p_end}
{phang2}{cmd:. margwals, asbalanced}{p_end}

{pstd}A more elaborated example on WALS average marginal effects{p_end}
{phang2}{cmd:. use "https://www.stata-press.com/data/r18/breathe", clear}{p_end}
{phang2}{cmd:. run "https://www.stata-press.com/data/r18/no2"}{p_end}
{phang2}{cmd:. display "$cc"}{p_end}
{phang2}{cmd:. display "$fc"}{p_end}
{phang2}{cmd:. vl create fc7 = fc - (grade)}{p_end}
{phang2}{cmd:. global ne "c.no2_class#i.grade c.no2_class#c.no2_class c.no2_class#c.no2_class#i.grade"}{p_end}
{phang2}{cmd:. wals react (c.no2_class i.grade) $ne $cc i.($fc7), rseed(12345)}{p_end}
{phang2}{cmd:. margwals, dydx(no2_class i.grade)}{p_end}
{phang2}{cmd:. margwals, dydx(no2_class) at(no2_class=(15(5)45))}{p_end}
{phang2}{cmd:. return list}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:margwals} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of observations{p_end}
{synopt:{cmd:r(k_predict)}}number of {opt predict()} options{p_end}
{synopt:{cmd:r(k_margins)}}number of terms in {it:marginlist}{p_end}
{synopt:{cmd:r(k_at)}}number of {opt at()} options{p_end}
{synopt:{cmd:r(level)}}confidence level of confidence intervals{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(cmd)}}{cmd:margwals}{p_end}
{synopt:{cmd:r(cmdline)}}command as typed{p_end}
{synopt:{cmd:r(est_cmd)}}{cmd:e(cmd)} from original estimation results{p_end}
{synopt:{cmd:r(est_cmdline)}}{cmd:e(cmdline)}
	from original estimation results{p_end}
{synopt:{cmd:r(title)}}title in output{p_end}
{synopt:{cmd:r(predict}{it:#}{cmd:_opts)}}the {it:#}th {cmd:predict()} option{p_end}
{synopt:{cmd:r(predict}{it:#}{cmd:_label)}}label from the {it:#}th {cmd:predict()} option{p_end}
{synopt:{cmd:r(expression)}}response expression{p_end}
{synopt:{cmd:r(xvars)}}{it:varlist} from {cmd:dydx()}, {cmd:dyex()},
					{cmd:eydx()}, or {cmd:eyex()}{p_end}
{synopt:{cmd:r(derivatives)}}"", "dy/dx", "dy/ex", "ey/dx", "ey/ex"{p_end}
{synopt:{cmd:r(over)}}{it:varlist} from {cmd:over()}{p_end}
{synopt:{cmd:r(within)}}{it:varlist} from {cmd:within()}{p_end}
{synopt:{cmd:r(atstats}{it:#}{cmd:)}}the {it:#}th {cmd:at()} specification
{p_end}
{synopt:{cmd:r(emptycells)}}{it:empspec} from {cmd:emptycells()}{p_end}

{p2col 5 20 24 2:Matrices}{p_end}
{synopt:{cmd:r(b)}}estimates{p_end}
{synopt:{cmd:r(_N)}}sample size corresponding to each margin estimate{p_end}
{synopt:{cmd:r(at)}}matrix of values from the {cmd:at()} options{p_end}
{synopt:{cmd:r(error)}}margin estimability codes;{break}
        {cmd:0} means estimable,{break}
        {cmd:8} means not estimable{p_end}
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

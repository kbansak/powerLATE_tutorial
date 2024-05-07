<!-- README.md is generated from README.Rmd. Please edit that file -->

# LATE Power Analysis Tutorial

------------------------------------------------------------------------

This page provides a tutorial on conducting power analysis for the Local
Average Treatment Effect (LATE)—also known as the Complier Average
Causal Effect (CACE)—and using the `powerLATE` package for that purpose.

------------------------------------------------------------------------

**Author:** [Kirk Bansak](https://www.kirkbansak.com/) (UC Berkeley)

**Date:** May 6, 2024

**Reference:** [Bansak, Kirk. “A generalized approach to power analysis
for local average treatment effects.” *Statistical Science* 35, no. 2
(2020).](https://doi.org/10.1214/19-STS732)

------------------------------------------------------------------------

# Contents

1.  Introduction

2.  Tutorial

# 1. Introduction

Performing a power analysis for the LATE can be an extremely important
step in one’s research design, just as it is when one’s interest is the
ATE (with perfect compliance). A power analysis can, for instance, allow
one to determine the necessary sample size or to understand the
magnitude of causal effects that can be detected given a fixed budget.
However, it turns out that analyzing power for estimators of the LATE is
much more complicated than for the ATE (or ITT). This may be surprising
given that the LATE can be viewed as a scaled version of the ITT—and in
some cases the variance of the most commonly used estimator of the LATE
(Wald IV estimator) is approximately equal to the variance of common ITT
estimators divided by the compliance rate. Indeed, in light of that
fact, one common recommendation is to use a heuristic power analysis
that focuses on the ITT even when one’s actual quantity of interest is
the LATE. Under some circumstances—i.e. given a suitable data-generating
process (DGP)—this heuristic approach can work well. However, it is
difficult and sometimes impossible to know in advance if one has such a
suitable DGP, while can result in the power analysis being misleading.

This problem is attributable to two factors. First, while the LATE is
indeed equivalent to the ITT scaled by compliance, the compliance rate
must also be estimated itself. Furthermore, the compliance rate and ITT
estimates will be correlated, and they can be correlated in any way
(positively or negatively) depending upon the unique features of the
DGP. This ends up having extremely important implications for the
statistical properties of the Wald IV estimator (or any estimator of the
LATE), including its variance and its power. As a result, as a number of
papers have shown (Jo 2002, Baiocchi et al. 2014, Bansak 2020), ITT and
related power analysis results diverge from the true power of tests of
the LATE, and sometimes dramatically so.

The second issue is that there is a proliferation of parameters that
factor into the variance (and hence power) of estimators of the LATE
compared to the ITT or ATE. For any power analysis, one can consider the
underlying parameters that factor into an estimator’s variance.
Considering the problem from a causal perspective can be done by
thinking about the relevant parameters in terms of distributions of
potential outcomes. For the ITT or ATE, there are two such “distribution
parameters”. In contrast, as shown by Bansak (2020), there are nine for
the LATE—related to the distributions of the treated and control
potential outcomes for the compliers, always-takers, and never-takers.
This severely complicates the power analysis enterprise and calls into
question standard approaches. It challenges the reliability of many
analytical instrumental-variable power analyses, which make simplifying
assumptions that can be unrealistic and/or require supplying parameters
one may not have good information on. Meanwhile, simulation-based
approaches require specifying all the parameters, or making modeling
choices that abstract some of them away; this leads to combinatorial
overload and/or potentially unreliable assumptions. The end result is
that, except in special high-information circumstances, results of the
power analysis can be misleading.

As a solution, Bansak (2020) proposes a bounds-based approach that does
not require specifying or making any assumptions about the distribution
parameters and provides conservative (though still practical) worst-case
power bounds. At a conceptual level, one can think of the method as
providing a single formula-based answer that captures the worst-case
scenario out of a countless number of simulations one would otherwise
need to run. In doing so, the method provides conservative answers for
the required sample size and minimum-detectable effects. The method
allows for the LATE to be specified either in terms of a raw effect or a
standardized effect size. The benefit of the former is interpretability,
but the cost is that it requires specifying one distributional
parameter—a similar case with standard ATE power analyses. That
distributional parameter can be dispensed with, however, if one uses the
standardized effect size. In addition, the method also allows for
further narrowing of the bounds (i.e. giving less conservative answers
and/or increasing the power) through the the introduction of an optional
assumption and/or the inclusion of covariates.

The method is implemented in the `R` package `powerLATE`. A tutorial is
provided below, with more information provided at
`https://github.com/kbansak/powerLATE`.

# 2. Tutorial

### 2.0. Getting started

A common (and perhaps the most common) reason researchers will perform a
power analysis is to determine the necessary sample size for their study
design to be able to reject the null hypothesis of no causal effect.
This tutorial will focus on that use case. Note, however, that the
method can also be used to calculate minimum detectable effects or
simply the power itself. Examples of these other use cases, which follow
a similar flow as the present use case, can be found in the Examples
section here: `https://github.com/kbansak/powerLATE`.

As mentioned above, the method can be implemented using the `powerLATE`
package in `R`.

``` r
install.packages("powerLATE")
```

``` r
library(powerLATE)
```

    ## powerLATE: Generalized Power Analysis for LATE
    ## Version: 0.1.1
    ## Reference: Bansak, K. (2020). A Generalized Approach to Power Analysis for Local Average Treatment Effects. Statistical Science, 35(2), 254-271.

As with any power analysis focused on sample size requirements, we must
specify a hypothetical magnitude/size of the treatment effect that we
are interested in. There are various ways a researcher can make this
decision. One possibility is to use the smallest effect that is deemed
to be of policy or theoretical interest. Another possibility is to use a
best guess of what the treatment effect will be based on prior studies
or evidence. Like with power analyses for the ATE, there are two ways in
which one can specify the magnitude of the treatment effect:

1.  The first is to specify the treatment effect in absolute terms. The
    LATE in absolute terms will be denoted here as *τ*. And as with a
    standard ATE power analysis, this also requires specifying the
    within-assignment group standard deviation of the outcome (which we
    denote here by *ω*). This requires some *a priori* empirical
    knowledge of the phenomenon of interest, and in practice researchers
    will commonly use the standard deviation of the outcome that been
    measured “in the wild”—i.e. in the absence of the treatment—based on
    either a pilot data collection or a related dataset. (More
    explicitly, the standard deviation of the outcome measured
    previously in the absence of the treatment corresponds to the
    outcome standard deviation for the control assignment group, which
    we take as an approximation of the pooled assignment group standard
    deviation.)

2.  Alternatively, instead of specifying the absolute magnitude of the
    LATE (*τ*) and within-group outcome standard deviation (*ω*), we can
    instead specify just the standardized LATE “effect size,” which we
    denote by *κ*. Popularized by psychologist/statistician Jacob Cohen
    in the 1980s (Cohen 1988), effect sizes refer to effects scaled by a
    reference standard deviation of the outcome. Here, as in standard
    ATE power analysis, the reference is the within-group standard
    deviation. The benefit of using effect sizes instead of absolute
    effects is the ability to avoid specifying any dispersion
    parameters, which one may not *a priori* have any reliable knowledge
    about. This does, however, require one to have an understanding of
    what constitutes a meaningful effect size based on theory, empirics,
    or precedent within one’s field.

As with a standard ATE power analysis, one also needs to specify the
assignment probability (i.e. proportion of units that will be assigned
to treatment in the study design, denoted here by *p*<sub>*z*</sub>),
the statistical power one desires (a conventional default is 0.8), and
the significance level *α* that one plans to use (a conventional default
is 0.05). Everything discussed so far that is needed for this LATE power
analysis is the same as what is needed for standard ATE power analyses.
Finally, however, there is one additional parameter needed that is
unique to the LATE: one must also specify an hypothesized compliance
rate (denoted by *π*). In practice, researchers will often specify a
range of possible values of *π* corresponding to different scenarios,
ideally informed by preliminary investigations or a pilot study.

Putting this all together, the LATE power analysis for computing a
sample size requires:

-   {*τ*, *ω*} OR *κ*
    -   This is the hypothesized effect along with the within-group
        outcome standard deviation OR just the hypothesized effect size.
    -   Same requirement for a standard ATE power analysis.
-   *p*<sub>*z*</sub>
    -   This is the proportion of units that will be assigned to the
        treatment in the study design.
    -   Same requirement for a standard ATE power analysis.
-   power (1 − *β*)
    -   This is the statistical power one desires.
    -   Same requirement for a standard ATE power analysis.
-   *α*
    -   This is the significance level (type I error tolerance) one
        plans to use.
    -   Same requirement for a standard ATE power analysis.
-   *π*
    -   This is the hypothesized compliance rate (or range of compliance
        rates).
    -   The one unique requirement here.

Below are examples using the `powerLATE` function, which presumes the
employment of the Wald IV estimator. We will use the conventional value
of 0.8 for the power (the `Power` argument in the function), and the
conventional value of 0.05 for the significance level (`sig.level`).
Further, in this example we will plan our study design such that the
probability of being assigned to the treatment (`pZ`) is 0.5. We will
specify a range of compliance values (`pi`) that we believe are all
plausible given subject matter expertise and/or external empirical
evidence—for illustration here it will include {0.6, 0.7, 0.8}.

### 2.1. Example 1: Calculating required sample size (using absolute effect magnitude)

We will begin by specifying an hypothesized LATE in absolute terms
(`tau`) along with a value of the within-group standard deviation
(`omega`). In this example, we have determined—based on subject matter
expertise, empirical knowledge about the outcome of interest, and
consultations with policymakers/experts—that there is a minimum effect
magnitude that is worthwhile. Anything smaller than that value is “too
small” to care, whereas if the effect is equal to or larger than that
value, we want to be sure we can detect it. This is the value at which
we set *τ*, and for illustration here we will set it to 25. (As a side
note, another way in which researchers will set *τ* is by performing
non-experimental analyses on existing data to try to approximate a
likely magnitude of the treatment effect.) To input a value for *ω*, we
can imagine that we have access to pre-study data in our context of
interest (or a similar context), and hence we can calculate the standard
deviation of the outcome for these units in the absence of the
treatment. For illustration here we will set it to 125.

We now have all the ingredients we need to compute our required sample
size(s). Finally, because we are specifying *τ* and *ω* (instead of just
*κ*), the function also requires that we set the `effect.size` argument
to `FALSE`. And now we run the function:

``` r
#Example 1
powerLATE(pi = c(0.6,0.7,0.8), power = 0.8, sig.level = 0.05,
          effect.size = FALSE, tau = 25, omega = 125)
```

    ## Power analysis for two-sided test that LATE equals zero
    ## 
    ##  pZ = 0.5
    ##  pi = Multiple values inputted (see table below)
    ##  tau = 25
    ##  omega = 125
    ##  Power = 0.8
    ##  sig.level  = 0.05 
    ## 
    ## Given these parameter values, the conservative (upper) bound for N (required sample size):
    ##   N        User-inputted pi
    ## 1 2543.037 0.6             
    ## 2 1838.766 0.7             
    ## 3 1377.969 0.8             
    ## 
    ## NOTE:  
    ## The Ordered-Means assumption is not being employed. If the user would like to make this assumption to narrow the bounds, set the argument assume.ord.means to TRUE.

As can be seen, since we provided three different compliance rate
possibilities, we are provided with three corresponding required sample
sizes. Note that the sample size pertains to the entire sample (i.e. not
the number per assignment group).

### 2.2. Example 2: Calculating required sample size (using standardized effect size)

Alternatively, we might also decide to specify our hypothesized LATE in
terms of a standardized effect size (*κ*). This is a particularly
attractive option if we do not have any reliable evidence about the
likely dispersion of the outcome (i.e. no good way to approximate *ω*).
To do so, we swap in the `kappa` argument while omitting the `tau` and
`omega` arguments, and set the `effect.size` argument to `TRUE`. For
instance, in this example, imagine that prior research on our phenomenon
of interest has argued that an effect size of 0.4 is economically
notable, and otherwise the treatment would not have sufficient policy
relevance. We will accordingly set `kappa` to 0.4 here:

``` r
#Example 2
powerLATE(pi = c(0.6,0.7,0.8), power = 0.8, sig.level = 0.05,
          effect.size = TRUE, kappa = 0.4)
```

    ## Power analysis for two-sided test that LATE equals zero
    ## 
    ##  pZ = 0.5
    ##  pi = Multiple values inputted (see table below)
    ##  kappa = 0.4
    ##  Power = 0.8
    ##  sig.level  = 0.05 
    ## 
    ## Given these parameter values, the conservative (upper) bound for N (required sample size):
    ##   N        User-inputted pi
    ## 1 733.4342 0.6             
    ## 2 523.0146 0.7             
    ## 3 384.5951 0.8             
    ## 
    ## NOTE:  
    ## The Ordered-Means assumption is not being employed. If the user would like to make this assumption to narrow the bounds, set the argument assume.ord.means to TRUE.

Again, we are provided with three required sample sizes that correspond
to the three possible compliance rates.

### 2.3. Example 3: Calculating minimum detectable effect size

As mentioned above, the `powerLATE` function can also be used to compute
minimum detectable effects (and effect sizes) or the power. Doing so
would follow the same flow as above, but would require inputting a value
for the sample size (`N`) in place of the other parameter(s). As one
brief example, the following will provide us with the minimum detectable
effect sizes if we have a fixed budget of 2000 units in our sample:

``` r
#Example 3
powerLATE(pi = c(0.6,0.7,0.8), power = 0.8, sig.level = 0.05,
          effect.size = TRUE, N = 2000)
```

    ## Power analysis for two-sided test that LATE equals zero
    ## 
    ##  pZ = 0.5
    ##  pi = Multiple values inputted (see table below)
    ##  N = 2000
    ##  Power = 0.8
    ##  sig.level  = 0.05 
    ## 
    ## Given these parameter values, the conservative (upper) bound for kappa (minimum detectable effect size):
    ##   kappa     User-inputted pi
    ## 1 0.2278494 0.6             
    ## 2 0.1912069 0.7             
    ## 3 0.1643345 0.8             
    ## 
    ## NOTE:  
    ## The Ordered-Means assumption is not being employed. If the user would like to make this assumption to narrow the bounds, set the argument assume.ord.means to TRUE.

Here, we are provided with three minimum detectable effect sizes that
correspond to the three possible compliance rates.

If we wanted to compute the minimum detectable effects (in absolute
terms), we would need to include `omega` as well. More on this can be
found in the Examples section here:
<https://github.com/kbansak/powerLATE>.

### 2.4 Example 4: Narrowing the bound

In the three examples shown so far, you may have noticed a note output
by the function stating, “The Ordered-Means assumption is not being
employed. If the user would like to make this assumption to narrow the
bounds, set the argument assume.ord.means to TRUE.” The Ordered-Means
assumption refers to an additional (optional) assumption that the user
can invoke, and in doing so make the results less conservative—for
instance, require smaller sample sizes. One must be careful invoking
this assumption: if it is wrong, and one mistakenly makes it, the power
analysis results will lead to an underpowered study. However, if one has
a strong theoretical argument or empirical evidence that the assumption
is met, then there are valuable practical gains from making the
assumption—such as saving time and/or money by requiring a smaller
sample. Before explaining what the assumption means substantively,
consider the following example, which shows how making the assumption
(setting `assume.ord.means = TRUE`) reduces the required sample size
relative to `Example 2` above:

``` r
#Example 4
powerLATE(pi = c(0.6,0.7,0.8), power = 0.8, sig.level = 0.05,
          effect.size = TRUE, kappa = 0.4, assume.ord.means = TRUE)
```

    ## Power analysis for two-sided test that LATE equals zero
    ## 
    ##  pZ = 0.5
    ##  pi = Multiple values inputted (see table below)
    ##  kappa = 0.4
    ##  Power = 0.8
    ##  sig.level  = 0.05 
    ## 
    ## Given these parameter values, the conservative (upper) bound for N (required sample size):
    ##   N        User-inputted pi
    ## 1 559.0147 0.6             
    ## 2 408.6223 0.7             
    ## 3 311.0119 0.8             
    ## 
    ## NOTE:  
    ## The Ordered-Means assumption is being employed. User should confirm that the assumption is reasonable in the context of interest.

As can be seen, sample sizes under the Ordered-Means assumption are
roughly 80% the size of the original sample sizes (i.e. a 20% savings in
sample size needs).

What does the Ordered-Means assumption mean, and when can it be reliably
invoked? Formally, the Ordered-Means assumption is that
*Ȳ*<sub>*N**T*</sub> ≤ *Ȳ*<sub>*C*</sub> ≤ *Ȳ*<sub>*A**T*</sub>, where
*Ȳ*<sub>*C*</sub>, *Ȳ*<sub>*N**T*</sub>, and *Ȳ*<sub>*A**T*</sub> denote
the expected realized outcome value for compliers, never-takers, and
always-takers. (And in the case of one-sided noncompliance, this can be
simplified to *Ȳ*<sub>*N**T*</sub> ≤ *Ȳ*<sub>*C*</sub> or
*Ȳ*<sub>*C*</sub> ≤ *Ȳ*<sub>*A**T*</sub>, depending upon the direction
of noncompliance.) As described in Bansak (2020):

> Roughly speaking, there are two factors to consider when assessing the
> plausibility of this assumption. The first relates to effect
> heterogeneity. Specifically, it should be the case that always-takers
> (never-takers) select into (out of) the treatment because treatment
> uptake for them is associated with effects that are larger (smaller)
> than the average treatment effect for the compliers, or at least
> similarly sized. For instance, in the case of a positive and
> beneficial treatment, we must expect the noncomplying study subjects
> to be sufficiently rational that they are selecting into (out of) the
> treatment in anticipation of a particularly good (bad) effect on their
> outcome. Alternatively, selection into and out of the treatment could
> also be made for arbitrary reasons that are uncorrelated with
> individual effects. The second factor relates to baseline outcome
> levels in the absence of the treatment. Specifically, we must expect
> that always-takers (never-takers) do not have baseline outcome levels
> that are particularly low (high) compared to that of the compliers.

More details, as well as ideas for how to make this assumption more
likely in the study design phase (thereby increasing the power of one’s
design), can be found in Bansak (2020).

### 2.5. Incorporating covariates and other examples

Finally, the `powerLATE` package also contains a function
`powerLATE.cov` that performs power analysis with covariate adjustment
(i.e. for the two-stage least squares estimator with covariate
adjustment). Implementation details can be found under the Examples
section here: `https://github.com/kbansak/powerLATE`

# References

1.  Baiocchi, Michael, Jing Cheng, and Dylan S. Small. “Instrumental
    variable methods for causal inference.” *Statistics in Medicine* 33,
    no. 13 (2014): 2297-2340.

2.  Bansak, Kirk. “A generalized approach to power analysis for local
    average treatment effects.” *Statistical Science* 35, no. 2 (2020):
    254-271.

3.  Cohen, Jacob. *Statistical Power Analysis for the Behavioral
    Sciences*. Hillsdale, NJ: Lawrence Erlbaum Associates (1988).

4.  Jo, Booil. “Statistical power in randomized intervention studies
    with noncompliance.” *Psychological Methods* 7, no. 2 (2002):
    178-193.

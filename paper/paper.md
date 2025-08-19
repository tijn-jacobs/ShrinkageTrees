---
title: 'ShrinkageTrees: an R package for Bayesian shrinkage trees'
tags:
  - R
  - C++
  - high-dimensional data
  - shrinkage priors
  - tree ensembles
  - causal inference
  - survival analysis
authors:
  - name: Tijn Jacobs
    orcid: 0009-0003-6188-9296
    affiliation: 1
affiliations:
 - name: Department of Mathematics, Vrije Universiteit Amsterdam, The Netherlands
   index: 1
   ror: 008xxew50
date: 31 August 2025
bibliography: paper.bib

---

# Summary

`ShrinkageTrees` provides Bayesian regression tree models with shrinkage priors
for high-dimensional prediction and causal inference. The package is tailored
for survival analysis with censored outcomes in settings where the number of
predictors exceeds the sample size ($p \gg n$). By shrinking irrelevant variables
instead of excluding them, it can estimate heterogeneous, non-linear treatment
effects while retaining important confounders.

In causal inference, valid estimation requires the *unconfoundedness*:
all covariates affecting both treatment assignment and the outcome must be
adjusted for. Methods that exclude covariates risk omitting important
confounders, which can bias treatment effect estimates. `ShrinkageTrees` protects
against such violations while still reducing noise in high-dimensional data by 
retaining all covariates and shrinking irrelevant ones toward zero.

`ShrinkageTrees` is relevant in fields such as genomics, epidemiology, and economics, 
where thousands of covariates may affect both treatment allocation
and outcomes. For example, in a genetic study of cancer patients, gene 
expression data can be used to adjust for confounding factors when estimating 
the effect of treatments such as radiation therapy on survival. An illustrative 
analysis of pancreatic cancer data is included and can be run using:

```r
demo("pdac_analysis", package = "ShrinkageTrees")
```


(DOUBLE:)
The package is aimed at applied researchers who analyse high-dimensional
datasets with censored outcomes, for example in genomics or clinical studies. It
is also useful for statisticians developing tree-based methodology who wish to
experiment with alternative priors on step heights. `ShrinkageTrees` is easy to
use in R [@R], integrates efficient C++ code via Rcpp [@Rcpp], and is available
on CRAN.


# Background

Let $T$ denote the non-negative survival time and $C$ the censoring time, with
observed follow-up $Y = \min(T,C)$ and censoring indicator $\delta \in \{0,1\}$.
Treatment assignment is indicated by $A \in \{0,1\}$, and covariates are denoted
by $X \in \mathbb{R}^p$. For prediction, `ShrinkageTrees` fits a single forest
model to capture the regression function $f(x)$ in
$\log T = f(x) + \varepsilon$, where $\varepsilon \sim \mathcal{N}(0,\sigma^2)$.
For causal inference, the outcome is decomposed into a prognostic component and
a treatment effect component,
$$
\log T(a) = f(x,\hat e(x)) + a \cdot \tau(x) + \varepsilon,
$$
where $\hat e(x)$ is the estimated probability of receiving treatment and 
$\tau(x)$ represents the
conditional treatment effect. Under standard assumptions (SUTVA,
unconfoundedness, positivity, and independent censoring), $\tau(x)$ identifies
the conditional average treatment effect (CATE),


For survival outcomes we work in the accelerated failure time (AFT) framework by
setting $Y = \log(T)$, and censored outcomes are handled through data
augmentation within the MCMC sampler.

The function $f$ is modelled by a Bayesian regression forest.
Each tree partitions the covariate space in, say $L$, subspaces.
to each partition, a step height is $h_\ell$ is assigned.
The step-heights are given a global--local shrinkage prior:



# Statement of need

For prediction, the package fits a single forest model to the outcome, enabling
accurate prediction of test observations or new data. For causal inference, it
implements a Bayesian Causal Forest [@bcf] decomposition of the outcome into a
prognostic component and a treatment-specific component, each modelled by its
own forest. This structure supports estimation of conditional average treatment
effects (CATEs) and allows flexible adjustment for confounders. For survival
outcomes, the package adopts the accelerated failure time (AFT) [@aft] framework, 
which models log-transformed survival times and naturally accommodates 
right-censoring.

The underlying methodology is described in [@Jacobs2025], where a horseshoe
prior is placed directly on the tree step heights to achieve adaptive
global–local shrinkage. This regularisation strategy preserves relevant signals
while reducing noise, improving performance in high-dimensional and sparse
settings. `ShrinkageTrees` is easy to use in R [@R], with efficient C++ 
integration via Rcpp [@Rcpp]. 



Tree-based methods are widely used because they flexibly capture non-linear 
effects and interactions. The Bayesian Additive Regression Trees (BART) [@BART] ...
Regularisation in these models is imposed through the tree structure, for example in Bayesian Causal Forests [@bcf], Dirichlet priors on splitting proportions [@dir], or sparse extensions [@sparsebcf]. In contrast,
`ShrinkageTrees` imposes shrinkage directly on the step heights through a
global–local prior. This provides adaptive regularisation while retaining all
covariates, offering better protection against violations of the unconfoundedness
assumption in causal inference. The package currently implements the horseshoe
prior [@horseshoe], with a general framework for scale mixture of normals priors.

For posterior computation, `ShrinkageTrees` uses an efficient reversible jump
Markov chain Monte Carlo sampler [@rjmcmc] with pseudo Gibbs proposals. To our
knowledge, there is no other general implementation in R that allows different
priors on step heights in BART-style models.

The package addresses prediction and causal inference with survival outcomes,
a setting where existing approaches are limited. For prediction, regularised Cox
models or survival SVMs focus on risk scores rather than flexible high-dimensional
non-linear effects. For causal inference, causal survival forests estimate
survival probabilities but are not adapted to high-dimensional settings. In
`ShrinkageTrees`, prediction is achieved by fitting a single forest model to the
outcome, while causal inference is based on a Bayesian Causal Forest decomposition
[@bcf] into prognostic and treatment-specific components. For survival data, the
package adopts the accelerated failure time (AFT) framework [@aft], which
log-transforms survival times and naturally accommodates right-censoring.

The methodology is detailed in [@Jacobs2025]. `ShrinkageTrees` targets applied
researchers in fields such as genomics, epidemiology, and biostatistics, where
datasets often contain thousands of predictors and censored outcomes. It is also
useful for statisticians developing tree-based methods who wish to experiment
with alternative priors on step heights. The package is easy to use in R [@R], 
with efficient C++ integration via Rcpp [@Rcpp], and is available on CRAN.

# Statement of need

Estimating treatment effects in high-dimensional data is challenging, especially
when outcomes are censored survival times. Many fields—such as genomics,
epidemiology, and biostatistics—routinely collect data with thousands of
covariates and relatively few observations. Valid causal inference requires the
unconfoundedness assumption: all covariates that affect both treatment assignment
and the outcome must be adjusted for. Methods that enforce sparsity by excluding
covariates risk omitting important confounders and producing biased estimates.
Existing approaches such as regularised Cox models or survival SVMs are designed
for prediction but not for high-dimensional causal inference, while causal
survival forests estimate survival probabilities but do not adapt well to
$p \gg n$ problems.

`ShrinkageTrees` addresses this gap by combining tree ensembles with global–local
shrinkage priors. For prediction, the package fits a single forest model to the
outcome. For causal inference, it implements a Bayesian Causal Forest
decomposition of the outcome into prognostic and treatment-specific components,
each modelled by its own forest. For survival data, it adopts the accelerated
failure time (AFT) framework, which models log-transformed survival times and
naturally accommodates right-censoring. Shrinkage is imposed directly on the
tree step heights through a horseshoe prior, reducing noise while retaining all
covariates to protect against violations of unconfoundedness. Posterior
inference is performed with an efficient reversible jump MCMC algorithm that
allows flexible extensions to other priors.




# Contributions

The C++ backend of `ShrinkageTrees` is designed to be modular, making it easy to
extend the package with other global–local shrinkage priors beyond the current
horseshoe implementations. Contributions are welcome: feel free to open an issue
to suggest new features or report bugs. If you wish to contribute code, you can
fork the repository, implement changes, and submit a pull request. Please ensure
that all tests pass, for example using `devtools::check()`, before submitting.

# Acknowledgements

Funded by the European Union. Views and opinions expressed are however those of 
the author(s) only and do not necessarily reflect those of the European Union or
the European Research Council Executive Agency. Neither the European Union nor
the granting authority can be held responsible for them. This work is supported 
by ERC grant BayCause, nr. 101074802.

# References
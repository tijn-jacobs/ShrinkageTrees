---
title: 'ShrinkageTrees: an R package for causal inference and prediction with Bayesian shrinkage trees in high-dimensions'
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

`ShrinkageTrees` provides Bayesian regression tree models with shrinkage
priors for high-dimensional prediction and causal inference. 
Unlike existing BART-based methods, it applies shrinkage directly to the step heights, offering flexible global–local priors. The package supports survival analysis with censored outcomes in settings where the
number of predictors exceeds the sample size ($p \gg n$). By shrinking
irrelevant variables instead of excluding them, it can estimate heterogeneous,
non-linear treatment effects while retaining important confounders.

In causal inference, valid estimation requires the *unconfoundedness*
assumption: all covariates affecting both treatment assignment and the outcome
must be accounted for. Methods that exclude covariates risk omitting important
confounders, which can bias treatment effect estimates. `ShrinkageTrees`
protects against such violations while still reducing noise in high-dimensional
data by retaining all covariates and shrinking irrelevant ones toward zero.

The package is aimed at applied researchers who analyse high-dimensional
datasets with censored outcomes, for example in genomics or clinical studies. 
For example, in a genetic study of cancer patients,
gene expression data can be used to adjust for confounding factors when
estimating the effect of treatments such as radiation therapy on survival. An
illustrative analysis of pancreatic cancer data is included and can be run
using:

```r
demo("pdac_analysis", package = "ShrinkageTrees")
```

`ShrinkageTrees` is easy to use in R [@R] and integrates efficient C++ code via Rcpp
[@Rcpp].
The methodology is described in [@Jacobs2025].


# Background

Let $T$ denote the---possibly censored---survival time.
Treatment assignment is indicated by $A \in \{0,1\}$, and covariates are denoted
by $x \in \mathbb{R}^p$. 
We specify an accelerated failure time model [@aft] for the survival times. 
The outcome is decomposed into a prognostic component $f$ and a treatment effect component $\tau$:
\begin{equation}
\log T = f(x) + A \cdot \tau(x) + \varepsilon,
\end{equation}
where $\varepsilon \sim \mathcal{N}(0,\sigma^2)$.
Under standard causal assumptions [@causalassumptions],
$\tau(x)$ identifies the *conditional average treatment effect*.
For prediction, `ShrinkageTrees` fits a single shrinkage
forest model to capture the regression function $f(x)$ in
$\log T = f(x) + \varepsilon$.
The censored outcomes are handled through data
augmentation within the Markov chain Monte Carlo sampler [@augmentation].


Both the prognostic function $f$ and the treatment effect function $\tau$ are
modelled by Bayesian shrinkage forests. Each tree partitions the covariate
space into $L$ regions. To each region $\ell\in\{1,2,...,L\}$ we assign a step height $h_\ell$,
which represents the contribution of that region to the overall prediction.
These step heights follow a global--local shrinkage prior:
$$
\begin{aligned}
h_\ell \mid \lambda_\ell^2, \gamma^2 &\sim \mathcal{N}(0, \lambda_\ell^2 \gamma^2), \\
\lambda_\ell^2 &\sim p(\lambda), \\
\gamma^2 &\sim p(\gamma),
\end{aligned}
$$
where the *global shrinkage parameter* $\gamma$ governs the overall degree
of shrinkage, encouraging all step heights to be small.
The *local shrinkage parameters* $\lambda_\ell$ allow individual step heights
$h_\ell$ to deviate from this global tendency when supported by the data. This
hierarchical structure suppresses noise while preserving meaningful signals in
both $f$ and $\tau$. 
The global--local prior encompasses a wide range of distributions [@globallocal].
In particular, choosing half-Cauchy priors for both $\gamma$ and $\lambda_\ell$ yields the horseshoe prior [@horseshoe].

# Statement of need


Estimating treatment effects in high-dimensional data is challenging, especially
when outcomes are censored survival times. Many fields---such as genomics,
epidemiology, and biostatistics---routinely collect data with thousands of
covariates and relatively few observations. Valid causal inference requires the
unconfoundedness assumption: all covariates that affect both treatment
assignment and the outcome must be adjusted for. Methods that enforce sparsity
by excluding covariates risk omitting important confounders and producing biased
estimates. 

Bayesian additive regression trees [@bart] and their ensemble extension are a popular tool for prediction and causal inference [@hill].
Regularisation in Bayesian regression trees is typically imposed through the tree structure. 
For example, regularisation can be imposed by assigning higher prior
probability to shallow trees [@bart], or by placing a Dirichlet prior on
the splitting probabilities. These strategies are implemented in, for instance,
the `BART` package.
In contrast, `ShrinkageTrees`
applies regularisation directly through the prior on the step heights.
`ShrinkageTrees` is, to our knowledge, the first software to allow flexible specification of shrinkage priors directly on the step heights.
This regularisation strategy preserves relevant signals
while reducing noise, improving performance in high-dimensional and sparse
settings.

Several extensions of BART for causal inference exist, such as Bayesian Causal Forests
[@bcf] implemented in the `bcf` package, and `SparseBCF`
[@sparsebcf] designed for high-dimensional covariates. However, these
methods are not designed for censored outcomes. The `AFTrees` package
[@Henderson] provides a BART implementation for survival data but fits a
single prognostic learner and does not address high-dimensional data. In contrast,
`ShrinkageTrees` combines high-dimensional causal inference with
survival outcomes, filling a gap in the current methodological and software
landscape. 

`ShrinkageTrees` extends the BART framework by introducing
flexible shrinkage priors on the step heights. The default implementation uses
the horseshoe prior, which provides adaptivity in high-dimensional settings.
This software fills a gap by accommodating censored outcomes, enabling both
prediction and causal inference for survival data with high-dimensional
covariates.



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
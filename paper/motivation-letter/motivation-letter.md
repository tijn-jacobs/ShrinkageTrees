---
output: pdf_document
fontsize: 12pt
---

\thispagestyle{empty}
\today

Editor  
The R Journal  
\bigskip

Dear Professor Tanaka,
\bigskip

Please consider our article "ShrinkageTrees: An R Package for Bayesian Tree Ensembles for Survival Analysis and Causal Inference" for publication in the R Journal.

The CRAN package `ShrinkageTrees` is aimed at applied researchers in biostatistics, clinical research, and epidemiology working with survival data. `BART` is currently the only CRAN option for Bayesian tree ensembles in this setting, and it does not support interval-censored outcomes, causal forests, or global--local shrinkage on the leaf parameters. `ShrinkageTrees` fills these gaps: it extends BART to right- and interval-censored outcomes via accelerated failure time models, provides the first R implementation of Bayesian causal forests for time-to-event endpoints, and introduces Horseshoe Forests, which place horseshoe priors on the leaf step heights for adaptive shrinkage in high-dimensional settings.

In line with the R Journal's software focus, the paper emphasises design and usage over methodology. We motivate the user-facing interface, S3 classes with `print`, `summary`, `predict`, and `plot` methods, an `Rcpp` sampler, integration with `ggplot2` and `coda`, and multi-chain parallel sampling, then illustrate a typical workflow on high-dimensional ovarian cancer survival data. The package is on CRAN, tested with `testthat`, fully documented with examples, accompanied by a vignette, and developed openly on GitHub. All code and data to reproduce the paper are included.

\bigskip
\bigskip

Regards,

Tijn Jacobs  
Department of Mathematics  
Vrije Universiteit Amsterdam  
Amsterdam, The Netherlands  
t.jacobs@vu.nl

\bigskip

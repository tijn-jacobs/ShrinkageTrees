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

Please consider our article titled "ShrinkageTrees: An R Package for Bayesian Tree Ensembles for Survival Analysis and Causal Inference" for publication in the R Journal.

Bayesian tree ensembles for survival data are in high demand across biostatistics, clinical research, and epidemiology, yet existing R packages cover only part of the landscape. The `BART` package is currently the only CRAN package that implements Bayesian tree ensembles for survival outcomes, and it does not support interval-censored data, causal forest models, or global--local shrinkage priors on the leaf parameters. Our package `ShrinkageTrees` addresses all three gaps. It extends BART to right-censored and interval-censored survival outcomes via accelerated failure time models, provides the first implementation of Bayesian causal forests for time-to-event endpoints in R, and implements Horseshoe Forests --- tree ensembles with horseshoe priors on the leaf step heights for adaptive shrinkage in high-dimensional settings. The package also supports Dirichlet splitting priors (DART) and allows both forms of regularisation to be combined within a single model.

Beyond the methodology, we have invested in making `ShrinkageTrees` a well-designed R package. It provides S3 classes with `print`, `summary`, `predict`, and `plot` methods, integrates with `ggplot2` for visualisation and `coda` for MCMC diagnostics, and supports multi-chain parallel sampling out of the box. The accompanying paper is dedicated to highlighting these software design choices and demonstrating the package through a worked example on high-dimensional ovarian cancer survival data.

We believe this article is a good fit for the R Journal because it introduces a package that fills a clear gap in the R ecosystem for Bayesian non-parametric survival analysis and causal inference, two areas of active methodological and applied interest.

\bigskip
\bigskip

Regards,

Tijn Jacobs  
Department of Mathematics  
Vrije Universiteit Amsterdam  
Amsterdam, The Netherlands  
t.jacobs@vu.nl

\bigskip

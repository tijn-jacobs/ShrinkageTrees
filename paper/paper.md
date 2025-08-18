---
title: 'ShrinkageTrees: an R package for Bayesian shrinkage trees'
tags:
  - R
  - high-dimensional data
  - shrinkage priors
  - tree ensembles
  - causal inference
  - survival analysis
authors:
  - name: Tijn Jacobs
    orcid: 0009-0003-6188-9296
    equal-contrib: true
    affiliation: 1
affiliations:
 - name: Department of Mathematics, Vrije Universiteit Amsterdam, The Netherlands
   index: 1
   ror: 008xxew50
date: 31 August 2025
bibliography: paper.bib

---

# Summary

`ShrinkageTrees` performs high-dimensional causal inference and prediction using
Bayesian regression trees with shrinkage priors. The package is specifically
tailored for survival analysis with censored data in high-dimensional settings,
where the number of predictors can exceed the sample size ($p>n$).

For prediction, the package fits a single forest model to the outcome, enabling
accurate prediction of test observations or new data. For causal inference, it
implements a Bayesian Causal Forest decomposition of the outcome into a
prognostic component and a treatment-specific component, each modelled by its
own forest. This structure supports estimation of conditional average treatment
effects (CATEs) and allows flexible adjustment for confounders. For survival
outcomes, the package adopts the accelerated failure time (AFT) framework, which
models log-transformed survival times and naturally accommodates 
right-censoring.

The underlying methodology is described in Jacobs (2025), where a horseshoe
prior is placed directly on the tree step heights to achieve adaptive
global–local shrinkage. This regularisation strategy preserves relevant signals
while reducing noise, improving performance in high-dimensional and sparse
settings. `ShrinkageTrees` is easy to use in R [@R], with efficient C++ 
integration via Rcpp [@Rcpp]. It also includes an illustrative analysis of 
pancreatic cancer data, which can be run via 
`demo("pdac_analysis", package = "ShrinkageTrees")`.

# Statement of need

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

Funded by the European Union. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them. This work is supported by ERC grant BayCause, nr. 101074802.

# References
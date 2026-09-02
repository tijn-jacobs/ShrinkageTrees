# Simulations

Simulation code and results for the two papers that accompany this package.
Each has its own folder, because the two use different data-generating
processes, different comparators, and different versions of the package.

## `r-journal-paper/`

Simulations for the software paper:

> Jacobs, T. *ShrinkageTrees: An R Package for Bayesian Tree Ensembles for
> Survival Analysis and Causal Inference.*
> [arXiv:2606.12317](https://arxiv.org/abs/2606.12317)

- `ovarian_analysis.R` -- the worked example (Sections 2 and 4.1-4.3), the
  manuscript's analysis chunks verbatim.
- `simulate_interval_censored.R` -- the interval-censored simulation study
  (Section 4.4): three priors, `p` in {50, 500, 5000}, 1000 replicates.
- `benchmark_timings.R` -- the wall-clock timing benchmark (Section 5.5).

All three require **ShrinkageTrees >= 2.1.0** and check the version on startup.
Earlier releases contain two defects in the horseshoe global update and use a
different default for `k`; see `NEWS.md`.

## `bayesian-analysis-paper/`

Simulations for the methodological paper:

> Jacobs, T., van Wieringen, W. N., & van der Pas, S. L. (2026). *Horseshoe
> Forests for High-Dimensional Causal Survival Analysis.* Bayesian Analysis,
> advance publication, 1-30. <https://doi.org/10.1214/26-BA1603>

See `bayesian-analysis-paper/README.md` for the folder layout. These scripts
predate version 2.1.0 and were run against the package as it stood at the time
of that submission.

# R Journal paper: simulations and worked example

Everything needed to regenerate the results in

> Jacobs, T. *ShrinkageTrees: An R Package for Bayesian Tree Ensembles for
> Survival Analysis and Causal Inference.*
> [arXiv:2606.12317](https://arxiv.org/abs/2606.12317)

All scripts require **ShrinkageTrees >= 2.1.0** and stop if an older version is
installed. Version 2.1.0 corrects two defects in the horseshoe global update and
recalibrates the default `k`; results from earlier releases do not reproduce.
See `NEWS.md` in the package root.

| file | produces |
|---|---|
| `ovarian_analysis.R` | the worked example, Sections 2 and 4.1-4.3 |
| `simulate_interval_censored.R` | the simulation study, Section 4.4 |
| `benchmark_timings.R` | the timing figure, Section 5.5 |

## `ovarian_analysis.R`

The worked example on the `ovarian` data: standard BART and a Horseshoe Forest
with out-of-sample prediction metrics, propensity scores, the causal forest, and
a posterior projection of the treatment-effect surface.

The analysis chunks of `ShrinkageTrees.Rmd` extracted verbatim and in the order
the manuscript runs them, so the numbers match. One deviation: the two tables
are rendered with `knitr::kable()` in the manuscript and printed as plain data
frames here.

```sh
Rscript ovarian_analysis.R
```

`ovarian` is semi-synthetic, so `ovarian_truth` carries the quantities the
models are trying to recover. RMSE, coverage and interval width are therefore
measured against the truth rather than against the observed outcome. The
generating script is `examples/semi-synthesise-ovarian.R` in the package root.

## `simulate_interval_censored.R`

The interval-censored simulation study, Section 4.4. Three priors
(`SurvivalBART()`, `SurvivalDART()`, `HorseTrees()`), `p` in {50, 500, 5000},
`n = 200`, 1000 replicates per setting, evaluated on training and test samples.

```sh
Rscript simulate_interval_censored.R              # full run
Rscript simulate_interval_censored.R --p 500      # one dimension only
Rscript simulate_interval_censored.R --cores 8
Rscript simulate_interval_censored.R --out DIR    # write elsewhere
Rscript simulate_interval_censored.R 192          # bare core count
```

Writes one file, `simulate_interval_censored_output.rds`, with all dimensions
combined. It goes to `$TMPDIR` when that is set, as cluster jobs do, and
otherwise to `outputs/` beside the script. `$TMPDIR` is wiped when the job ends,
so copy the file somewhere permanent before the job exits.

On a cluster the core count comes from the scheduler (`SLURM_CPUS_PER_TASK` and
friends) or from a bare positional argument. `OMP_NUM_THREADS` is deliberately
ignored: job scripts set it to 1 to stop BLAS threading inside each worker,
which says nothing about how many workers to run.

Replicate `i` at dimension `p` uses seed `2026 + 1000000 * index(p) + i`, set
inside the worker and keyed to the value of `p` rather than its position in the
loop. Results therefore do not depend on the number of cores, on the scheduling
order, or on whether the dimensions are run together or separately: `--p 500`
alone reproduces exactly the `p = 500` rows of a full run.

The full run is expensive: 9000 fits, each with 5000 burn-in and 5000 posterior
draws on 200 trees, so roughly 110 core-hours. That is about 40 minutes on 192
cores and a day on a laptop. Run one dimension first to check the pipeline.

## `benchmark_timings.R`

The wall-clock timing benchmark, Section 5.5. Times `SurvivalBART()`,
`SurvivalDART()`, `HorseTrees()` and `BART::mc.abart()` over `n` in
{100, 250, 500, 1000, 2000} at fixed `p = 100`, `m = 200`, 100 posterior draws
and `n_chains = 4`, with 100 replicates per point.

```sh
Rscript benchmark_timings.R              # full run, about an hour
Rscript benchmark_timings.R --smoke      # one replicate per n
Rscript benchmark_timings.R --replot     # redraw the figure, no re-timing
```

Writes `outputs/benchmark_timings.rds`, `outputs/benchmark_timings.csv` and
`outputs/benchmark_scaling.pdf`. The `.rds` carries `sessionInfo()` as an
attribute, so a rerun elsewhere can be compared against the run behind the
published figure.

`--replot` reads the stored `.rds` and rewrites the PDF without fitting
anything, which is what you want when only the figure needs another pass. It is
the one part of the script worth running twice, so it lives behind a flag
rather than in a separate file.

`BART::mc.abart()` splits `ndpost` across cores while ShrinkageTrees treats
`N_post` as per chain, so the script multiplies it up. Without that correction
`abart()` does a quarter of the posterior work and appears correspondingly
faster.

**Timings are hardware dependent.** The manuscript reports a 2024 MacBook Air
(Apple M3, 16 GB) running R 4.5.2, ShrinkageTrees 2.1.0 and BART 2.9.10.
Absolute seconds will differ elsewhere, and on a loaded or thermally throttled
machine they will differ substantially. The claim the figure supports is the
relative ordering and the approximately linear scaling in `n`, not the absolute
values. This is why the figure is generated here and included in the manuscript
as a static PDF rather than rebuilt during the knit: a paper whose numbers
change with the machine that compiles it is harder to check, not easier.

`BART` and `ggplot2` are optional. Without `BART` the `abart()` comparison is
skipped; without `ggplot2` the timings are still written and only the figure is
skipped.

# R Journal paper: simulations

Everything needed to regenerate the simulation results and the timing figure in

> Jacobs, T. *ShrinkageTrees: An R Package for Bayesian Tree Ensembles for
> Survival Analysis and Causal Inference.* R Journal (submitted).

Both scripts require **ShrinkageTrees >= 2.1.0** and stop if an older version is
installed. Version 2.1.0 corrects two defects in the horseshoe global update and
recalibrates the default `k`; results from earlier releases do not reproduce.
See `NEWS.md` in the package root.

## `simulate_interval_censored.R`

The interval-censored simulation study, Section 4.4 and Table 3. Three priors
(`SurvivalBART()`, `SurvivalDART()`, `HorseTrees()`), `p` in {50, 500, 5000},
`n = 200`, 1000 replicates per setting, evaluated on training and test samples.

```sh
Rscript simulate_interval_censored.R              # full run
Rscript simulate_interval_censored.R --smoke      # 2 replicates, tiny MCMC
Rscript simulate_interval_censored.R --p 500      # one dimension only
Rscript simulate_interval_censored.R --cores 8
```

Writes a single file, `outputs/simulate_interval_censored_output.rds`, with all
dimensions combined. On a cluster it goes to `$TMPDIR` instead unless `--out` is
given.

Replicate `i` at dimension `p` uses seed `20250101 + 1000 * index(p) + i`, set
inside the worker and keyed to the value of `p` rather than its position in the
loop. Results therefore do not depend on the number of cores, on the scheduling
order, or on whether the dimensions are run together or separately: `--p 500`
alone reproduces exactly the `p = 500` rows of a full run.

The full run is expensive. It is 9000 fits, each with 5000 burn-in and 5000
posterior draws on 200 trees; budget a cluster or a long weekend. Use `--smoke`
first to confirm the pipeline works end to end.

## `benchmark_timings.R`

The wall-clock timing benchmark, Section 6.5 and Figure 6. Times
`SurvivalBART()`, `SurvivalDART()`, `HorseTrees()` and `BART::mc.abart()` over
`n` in {100, 250, 500, 1000, 2000} at fixed `p = 100`, `m = 200`, 1000 posterior
draws and `n_chains = 4`, with three replicates per point.

```sh
Rscript benchmark_timings.R              # full run
Rscript benchmark_timings.R --smoke      # two n values, one replicate
```

Writes `outputs/benchmark_timings.rds`, `outputs/benchmark_timings.csv` and
`outputs/benchmark_scaling.pdf`. The `.rds` carries `sessionInfo()` as an
attribute, so a rerun elsewhere can be compared against the run behind the
published figure.

**Timings are hardware dependent.** The manuscript reports a 2024 MacBook Air
(Apple M3, 16 GB) running R 4.5.2, ShrinkageTrees 2.1.0 and BART 2.9.10.
Absolute seconds will differ on other machines, and on a loaded or thermally
throttled machine they will differ substantially. The claim the figure supports
is the relative ordering and the approximately linear scaling in `n`, not the
absolute values. This is why the figure is generated here and included in the
manuscript as a static PDF rather than rebuilt during the knit: a paper whose
numbers change with the machine that compiles it is harder to check, not easier.

`BART` and `ggplot2` are optional. Without `BART` the `abart()` comparison is
skipped; without `ggplot2` the timings are still written and only the figure is
skipped.

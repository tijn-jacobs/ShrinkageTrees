# Simulations

This folder contains all simulation code and results for the methodological paper on Shrinkage Trees.

## Folder structure

### `evaluation functions/`
Shared evaluation/helper functions used across multiple simulation scripts:
- `evaluation_functions.R` -- general evaluation functions
- `evaluation_functions_main.R` -- evaluation functions for the main simulation
- `evaluation_functions_new.R` -- evaluation functions for the deeper simulation

### `main sim/`
Main simulation study comparing methods under linear and nonlinear treatment effect scenarios (three dimensions setting):
- `revision_main_linear_low.R` -- linear treatment effect simulation
- `revision_main_nonlinear_low.R` -- nonlinear treatment effect simulation
- `revision_main_deeper.R` -- deeper tree simulation
- `revision_main_plotting.R` -- plotting script for the main simulation results

### `simple/`
Simple/illustrative simulation scenarios and supplementary simulation studies:
- `revision_simple.R` -- simple simulation
- `revision_simple_deeper.R` -- simple simulation with deeper trees (includes combined plotting of 3a/3b/3c)
- `revision_simple_deeper3a.R`, `revision_simple_deeper3b.R`, `revision_simple_deeper3c.R` -- sub-simulations for the deeper simple scenario
- `revision_continuous.R` -- continuous treatment simulation
- `revision_continuous_grf.R` -- continuous treatment simulation with GRF
- `revision_main_plotting.R` -- alternative main simulation plotting script

### `RIC/`
Simulation study for the RIC (Relative Influence Criterion):
- `ric.R` -- RIC simulation script
- `evaluation_functions_ric.R` -- evaluation functions specific to RIC

### `PDAC data analysis/`
Real data analysis using the PDAC (Pancreatic Ductal Adenocarcinoma) dataset:
- `pdac_analysis_revision.R` -- PDAC data analysis script

### `runtime/`
Runtime comparison between methods:
- `runtime.R` -- runtime benchmarking script

### `review/`
Simulations used for responding to reviewer comments (not in the main manuscript):

- **`acic/`** -- ACIC (Atlantic Causal Inference Conference) benchmark simulations at various sample sizes (500, 1000, 5000)
- **`sensitivity/`** -- Sensitivity analyses for the PDAC data analysis (varying M, Z, and their proportions)
- **`Data analysis simulations/`** -- Simulation studies based on the PDAC data-generating process
- **`inclusion/`** -- Posterior inclusion probability plots
- **`old simulations/`** -- Simulations from a previous version of the manuscript:
  - `Main sim/` -- original main simulation (linear and nonlinear)
  - `Deeper Friedman/` -- deeper Friedman simulation study (low/medium/high, LD/HD)
  - `Robust sim/` -- robustness simulations (correlated, dense, misspecification)
  - `continuous.R` -- original continuous treatment simulation

## Output files

Simulation results are saved as `.rds` files alongside their corresponding R scripts. Plots are saved as `.png` or `.pdf` files in the same directories.

Simulation scripts that run on a cluster save output to `$TMPDIR` (using `Sys.getenv("TMPDIR")`); the resulting `.rds` files were then copied into this folder for local analysis and plotting.


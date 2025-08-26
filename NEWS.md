# ShrinkageTrees 1.0.0

🎉 First CRAN release of **ShrinkageTrees**!

This package provides Bayesian regression tree models with shrinkage priors, supporting:
- Continuous outcomes
- Binary outcomes
- Right-censored survival data

It includes four core functions:
- `HorseTrees()`: fits a single regression tree with a standard Horseshoe prior.
- `ShrinkageTrees()`: fits a single tree with customizable shrinkage priors.
- `CausalHorseForest()`: fits a causal forest using the standard Horseshoe prior.
- `CausalShrinkageForest()`: fits a flexible causal forest with user-defined shrinkage priors and tuning options.

The `...Trees` functions use a single learner to estimate the outcome model directly. In contrast, the `Causal...Forest` variants fit separate models for the treated and control regression function. This enables estimation of conditional average treatment effects (CATEs).

# ShrinkageTrees 1.0.1

- Fixed bugs in the demo script  
- Added a `CONTRIBUTING.md` file with guidelines for contributors  
- Corrected minor typos in the source code (non-functional changes)  
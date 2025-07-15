# ShrinkageTrees <img src="https://img.shields.io/badge/R%3E%3D-4.2-blue" alt="R >= 4.2"> ![License: MIT](https://img.shields.io/badge/license-MIT-green) ![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange) <img src="sticker/ShrinkageTrees_hex.png" align="right" width="140"/>



This package provides functions for fitting Horseshoe Trees, Causal Horseshoe Forests, and their more general counterparts: Shrinkage Trees and Causal Shrinkage Forests.  

These models allow for global-local shrinkage priors on tree step heights, enabling adaptive modeling in high-dimensional settings.

The functions can be used for:

1) High-dimensional prediction  
2) High-dimensional causal inference 
3) Estimaton of heterogeneous (conditional average) treatment effects

Supported outcome types: continuous, binary, and **right-censored survival times**.

The mathematical background and theoretical foundation for these models are described in the manuscript "Horseshoe Forests for High-Dimensional Causal Survival Analysis" by T. Jacobs, W. van Wieringen, and S. van der Pas (2025).


## ✨ Features

- Horseshoe, forest-wide horseshoe, empirical Bayes Horseshpe, and half-Cauchy priors
- Flexible tree-based non-linear modeling of the ATE and CATE
- Supports survival data with right-censoring (accelerated failure time model)
- Efficient C++ backend via Rcpp

## 📦 Installation

You cannot install the released version of ShrinkageTrees from [CRAN](https://CRAN.R-project.org) yet:

```r
# install.packages("ShrinkageTrees")
```

You can install the development version from [GitHub](https://github.com/tijn-jacobs/ShrinkageTrees):

```r
# Install devtools if not already installed
install.packages("devtools")
devtools::install_github("tijn-jacobs/ShrinkageTrees")
```


## 🚀 Example

```r
library(ShrinkageTrees)

set.seed(42)
n <- 100
p <- 1000

# Generate covariates
X <- matrix(runif(n * p), ncol = p)
X_treat <- X_control <- X
treatment <- rbinom(n, 1, X[, 1])
tau <- 1 + X[, 2]/2 - X[, 3]/3 + X[, 4]/4

# Generate survival times (on log-scale)
true_time <- X[, 1] + treatment * tau + rnorm(n)
censor_time <- log(rexp(n, rate = 0.05))
follow_up <- pmin(true_time, censor_time)
status <- as.integer(true_time <= censor_time)

# Fit a standard Causal Horseshoe Forest (without propensity score adjustment)
fit_horseshoe <- CausalHorseForest(
  y = follow_up,
  status = status,
  X_train_control = X_control,
  X_train_treat = X_treat,
  treatment_indicator_train = treatment,
  outcome_type = "right-censored",
  timescale = "lof",
  number_of_trees = 200,
  k = 0.1,
  N_post = 5000,
  N_burn = 5000,
  store_posterior_sample = TRUE
)

# Posterior mean CATEs
CATE_horseshoe <- colMeans(fit_horseshoe$train_predictions_sample_treat)

# Posteriors of the ATE
post_ATE_horseshoe <- rowMeans(fit_horseshoe$train_predictions_sample_treat)

# Posterior mean ATE
ATE_horseshoe <- mean(post_ATE_horseshoe)

# Plot the posterior of the ATE
```
![Posterior ATE plot](man/figures/posterior_ate_plot.png)


## 📄 Documentation

- In R: `?ShrinkageTrees`, `?HorseTrees`, `?CausalHorseForest`, and `?CausalShrinkageForest` for detailed help.
- Examples and parameter descriptions can be found in each function’s documentation.


## 🤝 Contributing

Contributions are welcome! Feel free to open an [issue](https://github.com/tijn-jacobs/ShrinkageTrees/issues) or submit a pull request. 
The software is designed to be flexible and modular, allowing for a wide variety of global-local shrinkage priors to be easily implemented and extended in future versions.



## 📄 License

This package is licensed under the [MIT License](https://cran.r-project.org/web/licenses/MIT).


## 💬 Acknowledgments

This project has received funding from the European Research Council (ERC) under the European Union’s Horizon Europe program under Grant agreement No. 101074082. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them

# ShrinkageTrees <img src="https://img.shields.io/badge/R%3E%3D-4.2-blue" alt="R >= 4.2"> ![License: MIT](https://img.shields.io/badge/license-MIT-green) ![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange)

This package provides functions for fitting Horseshoe Trees, Causal Horseshoe Forests, and their more general counterparts: Shrinkage Trees and Causal Shrinkage Forests.  

These models allow for **flexible global-local shrinkage priors** on tree step heights, enabling robust modeling in high-dimensional settings.

The functions can be used for:

✅ High-dimensional prediction  
✅ Causal inference of heterogeneous treatment effects given high-dimensional covariates  

Supported outcome types: **continuous**, **binary**, and **right-censored survival times**.


## ✨ Features

- Horseshoe, forest-wide horseshoe, empirical Bayes, and half-Cauchy priors
- Flexible tree and forest-based non-linear modeling
- Separate control and treatment trees for causal effect decomposition
- Supports survival data with right-censoring (accelerated failure time interpretation)
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

# Generate data
n <- 50
p <- 3
X <- matrix(runif(n * p), ncol = p)
y <- X[, 1] + rnorm(n)

# Fit a Shrinkage Tree model with standard horseshoe prior
fit_hs <- ShrinkageTrees(
  y = y,
  X_train = X,
  outcome_type = "continuous",
  number_of_trees = 10,
  prior_type = "horseshoe",
  local_hp = 0.1 / sqrt(10),
  global_hp = 0.1 / sqrt(10),
  N_post = 10,
  N_burn = 5,
  verbose = FALSE,
  seed = 1
)

# Fit with half-Cauchy prior
fit_hc <- ShrinkageTrees(
  y = y,
  X_train = X,
  outcome_type = "continuous",
  number_of_trees = 10,
  prior_type = "half-cauchy",
  local_hp = 1 / sqrt(10),
  N_post = 10,
  N_burn = 5,
  verbose = FALSE,
  seed = 1
)

# Compare posterior mean predictions
plot(fit_hs$train_predictions, type = "l", col = "steelblue", ylim = range(c(fit_hs$train_predictions, fit_hc$train_predictions)), ylab = "Prediction", main = "ShrinkageTrees predictions")
lines(fit_hc$train_predictions, col = "orange3")
legend("topright", legend = c("Horseshoe", "Half-Cauchy"), col = c("steelblue", "orange3"), lty = 1)
```


## 📄 Documentation

- In R: `?ShrinkageTrees`, `?HorseTrees`, `?CausalHorseForest`, and `?CausalShrinkageForest` for detailed help.
- Examples and parameter descriptions can be found in each function’s documentation.


## 🤝 Contributing

Contributions are welcome! Feel free to open an [issue](https://github.com/tijn-jacobs/ShrinkageTrees/issues) or submit a pull request. 
The software is designed to be flexible and modular, allowing for a wide variety of global-local shrinkage priors to be easily implemented and extended in future versions.



## 📄 License

This package is licensed under the **MIT License**.






## 💬 Acknowledgments

This project has received funding from the European Research Council (ERC) under the European Union’s Horizon Europe program under Grant agreement No. 101074082. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them

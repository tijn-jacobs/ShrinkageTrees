# ── Example: MATE vs. PATE (Bayesian bootstrap) ───────────────────────────────
#
# Fits a CausalHorseForest on a small simulated dataset with heterogeneous
# treatment effects, then compares the two ATE posteriors:
#
#   - MATE: mixed ATE, averaging CATEs with equal 1/n weights at each
#     iteration (conditions on the observed covariates).
#   - PATE: population ATE, reweighting CATEs with Dirichlet(1, ..., 1)
#     weights (integrates over F_X -- correct credible interval).
#
# The two posteriors are printed side by side and plotted on a shared axis
# so the CI widening from the Bayesian bootstrap is visible at a glance.
#
# Run from the package root:
#   Rscript examples/test-bayesian-bootstrap-ate.R

library(ShrinkageTrees)
library(ggplot2)

# ── 1. Simulated data with heterogeneous tau(x) ───────────────────────────────

set.seed(1)
n <- 300; p <- 5
X <- matrix(rnorm(n * p), ncol = p)
treatment <- rbinom(n, 1, 0.5)

mu_x  <- X[, 1] + 0.5 * X[, 2]
tau_x <- 1 + X[, 3] + 0.5 * X[, 4]          # heterogeneous; true E[tau] ≈ 1
y     <- mu_x + (treatment - 0.5) * tau_x + rnorm(n)

# ── 2. Fit ────────────────────────────────────────────────────────────────────

fit <- CausalHorseForest(
  y                         = y,
  X_train_control           = X,
  X_train_treat             = X,
  treatment_indicator_train = treatment,
  outcome_type              = "continuous",
  N_post                    = 1000,
  N_burn                    = 500,
  store_posterior_sample    = TRUE,
  verbose                   = FALSE
)

# ── 3. Side-by-side ATE posteriors ────────────────────────────────────────────

bb <- bayesian_bootstrap_ate(fit)

cat("\n== ATE posteriors ==\n")
cat(sprintf("  MATE:  %.3f   95%% CI: [%.3f, %.3f]   width = %.3f\n",
            bb$mate_mean, bb$mate_ci$lower, bb$mate_ci$upper,
            bb$mate_ci$upper - bb$mate_ci$lower))
cat(sprintf("  PATE:  %.3f   95%% CI: [%.3f, %.3f]   width = %.3f\n",
            bb$pate_mean, bb$pate_ci$lower, bb$pate_ci$upper,
            bb$pate_ci$upper - bb$pate_ci$lower))

# Same info via summary() with the bayesian_bootstrap toggle:
summary(fit, bayesian_bootstrap = TRUE)   # default — PATE
summary(fit, bayesian_bootstrap = FALSE)  # MATE

# ── 4. Overlay plot ───────────────────────────────────────────────────────────

draws <- rbind(
  data.frame(ate = bb$mate_samples, kind = "MATE"),
  data.frame(ate = bb$pate_samples, kind = "PATE")
)

ci_df <- rbind(
  data.frame(kind = "MATE",
             lower = bb$mate_ci$lower, upper = bb$mate_ci$upper),
  data.frame(kind = "PATE",
             lower = bb$pate_ci$lower, upper = bb$pate_ci$upper)
)

p_overlay <- ggplot(draws, aes(x = ate, fill = kind, colour = kind)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = ci_df, aes(xintercept = lower, colour = kind),
             linetype = "dashed") +
  geom_vline(data = ci_df, aes(xintercept = upper, colour = kind),
             linetype = "dashed") +
  scale_fill_manual(values = c(MATE = "grey40", PATE = "steelblue")) +
  scale_colour_manual(values = c(MATE = "grey40", PATE = "steelblue")) +
  labs(x = "ATE", y = "Posterior density",
       title = "MATE vs. PATE (Bayesian bootstrap) -- 95% CIs dashed",
       fill = NULL, colour = NULL) +
  theme_bw()

print(p_overlay)

# ── 5. Same via plot.CausalShrinkageForest() ──────────────────────────────────
# Identical effect via the built-in method, one plot per variant.

plot(fit, type = "ate", bayesian_bootstrap = TRUE)   # default
plot(fit, type = "ate", bayesian_bootstrap = FALSE)

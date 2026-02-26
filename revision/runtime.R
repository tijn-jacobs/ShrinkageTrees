## ============================================================
## Runtime comparison: BART vs Causal Horseshoe Forest (CHF)
## Local MacBook version with silent BART + progress updates
## ============================================================

suppressPackageStartupMessages({
  library(BART)
  library(ShrinkageTrees)
  library(foreach)
  library(doParallel)
  library(ggplot2)
})

## ------------------------------------------------------------
## CONFIGURATION
## ------------------------------------------------------------
n_cores <- 6          # <-- CHANGE THIS (e.g. 5 or 6)
R       <- 200        # repetitions per (DGP, m)
m_values <- seq(10, 250, by = 10)
ndpost  <- 1000
nskip   <- 1000

## ------------------------------------------------------------
## Parallel backend (used ONLY to silence BART output)
## ------------------------------------------------------------
cl <- makeCluster(n_cores)
registerDoParallel(cl)

cat("Runtime comparison: BART vs CHF\n")
cat("Cores used:", n_cores, "\n")
cat("Repetitions per setting:", R, "\n\n")

## ------------------------------------------------------------
## Wall-clock timing
## ------------------------------------------------------------
time_it <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - t0
}

## ------------------------------------------------------------
## Data-generating processes (5 DGPs)
## ------------------------------------------------------------
generate_data <- function(dgp_id, seed) {
  set.seed(seed)
  
  if (dgp_id == 1) {
    n <- 50;  p <- 5
  } else if (dgp_id == 2) {
    n <- 50;  p <- 50
  } else if (dgp_id == 3) {
    n <- 100; p <- 5
  } else if (dgp_id == 4) {
    n <- 100; p <- 100
  } else if (dgp_id == 5) {
    n <- 50;  p <- 250
  } else {
    stop("dgp_id must be 1..5")
  }
  
  X <- matrix(runif(n * p), n, p)
  a <- rbinom(n, 1, 0.5)
  k <- min(5L, p)
  y <- rowSums(X[, 1:k, drop = FALSE]) + a + rnorm(n)
  
  list(X = X, y = y, a = a, n = n, p = p)
}

## ------------------------------------------------------------
## Single runtime evaluation
## ------------------------------------------------------------
run_single_runtime <- function(method, X, y, a, m, ndpost, nskip) {
  
  if (method == "CHF") {
    return(
      time_it({
        ShrinkageTrees::CausalShrinkageForest(
          y = y,
          status = NULL,
          outcome_type = "continuous",
          X_train_control = X,
          X_train_treat   = X,
          treatment_indicator_train = a,
          number_of_trees_control = m,
          number_of_trees_treat   = m,
          prior_type_control = "horseshoe",
          prior_type_treat   = "horseshoe",
          local_hp_control  = 0.1,
          local_hp_treat    = 0.1,
          global_hp_control = 0.1,
          global_hp_treat   = 0.1,
          N_post = ndpost,
          N_burn = nskip,
          store_posterior_sample = FALSE,
          verbose = FALSE
        )
      })
    )
  }
  
  if (method == "BART") {
    return(
      time_it({
        sink(tempfile())     # silence BART
        on.exit(sink(), add = TRUE)
        
        BART::wbart(
          x.train = cbind(X, Z = a),
          y.train = y,
          ntree   = 2*m,
          ndpost  = ndpost,
          nskip   = nskip
        )
      })
    )
  }
  
  stop("Unknown method")
}

## ------------------------------------------------------------
## Task grid
## ------------------------------------------------------------
grid <- expand.grid(
  dgp_id = 1:5,
  m      = m_values,
  rep    = seq_len(R)
)

n_tasks <- nrow(grid)
chunk_size <- 25L
n_chunks <- ceiling(n_tasks / chunk_size)

cat("Total tasks:", n_tasks, "\n")
cat("Chunk size:", chunk_size, "\n\n")

## ------------------------------------------------------------
## Run simulation with progress updates (MASTER ONLY)
## ------------------------------------------------------------
runtime_results <- vector("list", n_chunks)

for (chunk in seq_len(n_chunks)) {
  
  idx_start <- (chunk - 1L) * chunk_size + 1L
  idx_end   <- min(chunk * chunk_size, n_tasks)
  idx_range <- idx_start:idx_end
  
  cat(sprintf(
    "[%s] Running tasks %d–%d of %d\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    idx_start, idx_end, n_tasks
  ))
  flush.console()
  
  runtime_results[[chunk]] <- foreach(
    i = idx_range,
    .combine  = rbind,
    .packages = c("BART", "ShrinkageTrees")
  ) %dopar% {
    
    g <- grid[i, ]
    
    dat <- generate_data(
      dgp_id = g$dgp_id,
      seed   = 10000L * g$dgp_id + g$rep
    )
    
    rbind(
      data.frame(
        dgp = g$dgp_id,
        n   = dat$n,
        p   = dat$p,
        m   = g$m,
        rep = g$rep,
        method = "CHF",
        runtime_sec = run_single_runtime(
          "CHF", dat$X, dat$y, dat$a, g$m, ndpost, nskip
        )
      ),
      data.frame(
        dgp = g$dgp_id,
        n   = dat$n,
        p   = dat$p,
        m   = g$m,
        rep = g$rep,
        method = "BART",
        runtime_sec = run_single_runtime(
          "BART", dat$X, dat$y, dat$a, g$m, ndpost, nskip
        )
      )
    )
  }
}

runtime_results <- do.call(rbind, runtime_results)

## ------------------------------------------------------------
## Aggregate over DGPs and repetitions
## ------------------------------------------------------------
runtime_summary <- aggregate(
  runtime_sec ~ m + method,
  data = runtime_results,
  FUN = function(x) c(mean = mean(x), sd = sd(x))
)

runtime_summary$runtime_mean <- runtime_summary$runtime_sec[, "mean"]
runtime_summary$runtime_sd   <- runtime_summary$runtime_sec[, "sd"]

## ------------------------------------------------------------
## Plot: runtime vs m
## ------------------------------------------------------------
p <- ggplot(runtime_summary,
            aes(x = m, y = runtime_mean,
                color = method, group = method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  geom_errorbar(
    aes(ymin = runtime_mean,
        ymax = runtime_mean),
    width = 8,
    alpha = 0.35
  ) +
  scale_color_manual(
    values = c("BART" = "#1f78b4",
               "CHF"  = "#ff3b1f")
  ) +
  labs(
    x = "Trees per forest (m)",
    y = "Average runtime (seconds)",
    color = "Method",
    title = "Runtime comparison: BART vs Causal Horseshoe Forest"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

print(p)

## ------------------------------------------------------------
## Save results
## ------------------------------------------------------------
saveRDS(
  list(raw = runtime_results, summary = runtime_summary),
  file = "runtime_comparison_BART_vs_CHF.rds"
)

cat("\nResults saved to runtime_comparison_BART_vs_CHF.rds\n")

## ------------------------------------------------------------
## Cleanup
## ------------------------------------------------------------
stopCluster(cl)
registerDoSEQ()








# Load the saved results
results <- readRDS("runtime_comparison_BART_vs_CHF.rds")

# Extract the individual components
raw_data <- results$raw
summary_data <- results$summary






# IQR cleaning because laptop was closed during runtime, and time kept being recorded

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# 1. Load the existing results
if(!file.exists("runtime_comparison_BART_vs_CHF.rds")) {
  stop("File not found! Please ensure the .rds file is in your working directory.")
}
results <- readRDS("runtime_comparison_BART_vs_CHF.rds")
raw_data <- results$raw

summary_iqr <- raw_data %>%
  group_by(m, method) %>%
  filter(runtime_sec < (quantile(runtime_sec, 0.75) + 1.5 * IQR(runtime_sec))) %>%
  summarise(
    runtime_mean = mean(runtime_sec),
    runtime_sd = sd(runtime_sec),
    .groups = "drop"
  ) %>%
  mutate(clean_method = "IQR Statistical Cleaning")


library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# Prepare data
# ------------------------------------------------------------
plot_data <- summary_iqr %>%
  mutate(
    method = factor(
      method,
      levels = c("BART", "CHF"),
      labels = c("BART", "Causal Horseshoe Forest")
    )
  )

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
p <- ggplot(plot_data,
            aes(x = m, y = runtime_mean,
                color = method, fill = method)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.6) +
  geom_ribbon(
    aes(
      ymin = runtime_mean - runtime_sd,
      ymax = runtime_mean + runtime_sd
    ),
    alpha = 0.20,
    linewidth = 0
  ) +
  scale_color_manual(
    values = c(
      "BART" = "#1f78b4",
      "Causal Horseshoe Forest" = "#e31a1c"
    )
  ) +
  scale_fill_manual(
    values = c(
      "BART" = "#1f78b4",
      "Causal Horseshoe Forest" = "#e31a1c"
    )
  ) +
  labs(
    x = "Trees per forest (m)",
    y = "Average runtime (seconds)",
    color = "Method",
    fill  = "Method",
    title = "Runtime comparison",
    subtitle = "Mean ± 1 SD across DGPs and repetitions"
  ) +
  theme_minimal(base_size = 15, base_family = "Times New Roman") +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )
print(p)

# ------------------------------------------------------------
# Export to PNG
# ------------------------------------------------------------
ggsave(
  filename = "runtime_comparison_BART_vs_CHF.png",
  plot = p,
  width = 9,
  height = 5,
  dpi = 300
)

cat("Plot saved as runtime_comparison_BART_vs_CHF.png\n")



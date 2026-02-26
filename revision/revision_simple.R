library(ShrinkageTrees)
library(BART)
library(SoftBart)
library(doParallel)
library(foreach)

data_gen <- function(
    n = 200,
    p = 100,
    sigma = 1,
    sparse = FALSE,
    sparsity = 0.05
) {
  
  # 1. Covariates
  X <- matrix(runif(n * p), nrow = n, ncol = p)
  
  # 2. Nonlinear Friedman signal (low-dimensional)
  f_x <- 10 * sin(pi * X[, 1] * X[, 2]) +
    20 * (X[, 3] - 0.5)^2 +
    10 * X[, 4] +
    5 * X[, 5]
  
  # 3. Optional sparse high-dimensional linear component
  if (sparse) {
    beta <- rbinom(p, 1, sparsity) * rnorm(p)
    f_x <- f_x + X %*% beta
  }
  
  # 4. Normalize signal
  f_x <- as.numeric(f_x)
  sd_Y <- sd(f_x)
  f_x <- f_x / sd_Y
  
  # 5. Add noise
  Y <- f_x + rnorm(n, 0, sigma)
  
  list(
    X = X,
    Y = as.numeric(Y),
    f_true = f_x
  )
}

compute_metrics <- function(f_hat, f_post, f_true) {
  
  rmse <- sqrt(mean((f_hat - f_true)^2))
  
  ci_lower <- apply(f_post, 2, quantile, 0.025)
  ci_upper <- apply(f_post, 2, quantile, 0.975)
  
  coverage <- mean(f_true >= ci_lower &
                     f_true <= ci_upper)
  
  length <- mean(ci_upper - ci_lower)
  
  c(RMSE = rmse,
    Coverage = coverage,
    Length = length)
}

fit_chf <- function(
    data,
    k = 0.1,
    N_post = 2500,
    N_burn = 2500,
    m = 200
) {
  
  # Optional: scale outcome (keep if your DGP scales signal)
  sd_Y <- 1
  Y_train <- data$Y / sd_Y
  
  st_fit <- ShrinkageTrees::ShrinkageTrees(
    y = Y_train,
    X_train = data$X,
    outcome_type = "continuous",
    number_of_trees = m,
    
    prior_type = "horseshoe",
    
    local_hp  = k / sqrt(m),
    global_hp = k / sqrt(m),
    
    N_post = N_post,
    N_burn = N_burn,
    
    store_posterior_sample = TRUE,
    verbose = FALSE
  )
  
  list(
    f_hat  = as.numeric(st_fit$train_predictions) * sd_Y,
    f_post = st_fit$train_predictions_sample * sd_Y,
    sigma_post = st_fit$sigma
  )
}


fit_bart <- function(
    data,
    N_post = 2500,
    N_burn = 2500
) {
  
  bart_fit <- capture.output(
    fit <- BART::wbart(
      x.train = data$X,
      y.train = data$Y,
      nskip   = N_burn,
      ndpost  = N_post,
      printevery = 999999999
    ),
    type = "output"
  )
  
  bart_fit <- fit
  
  # Posterior draws: [ndpost × n]
  f_post <- bart_fit$yhat.train
  
  f_hat <- colMeans(f_post)
  
  list(
    f_hat  = as.numeric(f_hat),
    f_post = f_post,
    sigma_post = bart_fit$sigma
  )
}

fit_dart <- function(
    data,
    N_post = 2500,
    N_burn = 2500
) {
  
  dart_fit <- capture.output(
    fit <- BART::wbart(
      x.train = data$X,
      y.train = data$Y,
      nskip   = N_burn,
      ndpost  = N_post,
      sparse = TRUE,
      printevery = 999999999
    ),
    type = "output"
  )
  
  dart_fit <- fit
  
  # Posterior draws: [ndpost × n]
  f_post <- dart_fit$yhat.train
  
  f_hat <- colMeans(f_post)
  
  list(
    f_hat  = as.numeric(f_hat),
    f_post = f_post,
    sigma_post = dart_fit$sigma
  )
}



fit_softbart <- function(
    data,
    m = 200,
    N_post = 2500,
    N_burn = 2500
) {
  
  # Hyperparameters
  hypers <- SoftBart::Hypers(
    X = data$X,
    Y = data$Y,
    num_tree = m
  )
  
  opts <- SoftBart::Opts(
    num_burn = N_burn,
    num_save = N_post
  )
  
  sb_fit <- SoftBart::softbart(
    X = data$X,
    Y = data$Y,
    X_test = data$X,   # predict on training set
    hypers = hypers,
    opts = opts,
    verbose = FALSE
  )
  
  # Posterior draws: [num_save × n]
  f_post <- sb_fit$y_hat_test
  f_hat  <- colMeans(f_post)
  
  list(
    f_hat  = as.numeric(f_hat),
    f_post = f_post,
    sigma_post = sb_fit$sigma
  )
}


run_one_sim <- function(
    seed,
    n = 100,
    p = 1000,
    sigma = 1,
    N_post = 2000,
    N_burn = 2000
) {
  
  set.seed(seed)
  
  dat <- data_gen(
    n = n,
    p = p,
    sigma = sigma,
    sparse = TRUE
  )
  
  # Fit models
  chf0.1_res   <- fit_chf(dat, k = 0.1, N_post = N_post, N_burn = N_burn)
  # chf0.2_res   <- fit_chf(dat, k = 0.2, N_post = N_post, N_burn = N_burn)
  # chf0.3_res   <- fit_chf(dat, k = 0.3, N_post = N_post, N_burn = N_burn)
  bart_res     <- fit_bart(dat, N_post = N_post, N_burn = N_burn)
  dart_res     <- fit_dart(dat, N_post = N_post, N_burn = N_burn)
  softbart_res <- fit_softbart(dat, N_post = N_post, N_burn = N_burn)
  
  # Metrics
  chf0.1_metrics      <- compute_metrics(chf0.1_res$f_hat,
                                         chf0.1_res$f_post,
                                         dat$f_true)
  
  #   chf0.2_metrics      <- compute_metrics(chf0.2_res$f_hat,
  #                                          chf0.2_res$f_post,
  #                                          dat$f_true)
  # 
  #   chf0.3_metrics      <- compute_metrics(chf0.3_res$f_hat,
  #                                          chf0.3_res$f_post,
  #                                          dat$f_true)
  
  bart_metrics     <- compute_metrics(bart_res$f_hat,
                                      bart_res$f_post,
                                      dat$f_true)
  
  dart_metrics     <- compute_metrics(dart_res$f_hat,
                                      dart_res$f_post,
                                      dat$f_true)
  
  softbart_metrics <- compute_metrics(softbart_res$f_hat,
                                      softbart_res$f_post,
                                      dat$f_true)
  
  rbind(
    CHF0.1       = chf0.1_metrics,
    # CHF0.2       = chf0.2_metrics,
    # CHF0.3       = chf0.3_metrics,
    BART      = bart_metrics,
    DART      = dart_metrics,
    SoftBART  = softbart_metrics
  )
}


run_simulation_parallel <- function(
    n_rep = 50,
    n = 100,
    p = 1000,
    sigma = 1,
    N_post = 2000,
    N_burn = 2000,
    n_cores = parallel::detectCores() - 1
) {
  
  results <- foreach(
    i = 1:n_rep,
    .combine = 'c',
    .packages = c("ShrinkageTrees",
                  "BART",
                  "stochtree",
                  "SoftBart"),
    .export = c(
      "run_one_sim",
      "data_gen",
      "compute_metrics",
      "fit_chf",
      "fit_bart",
      "fit_dart",
      "fit_softbart"
    )
  ) %dopar% {
    
    run_one_sim(
      seed = i,
      n = n,
      p = p,
      sigma = sigma,
      N_post = N_post,
      N_burn = N_burn
    )
  }
  
  results
}






# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  num_cores <- as.integer(args[1]) - 1
} else {
  num_cores <- parallel::detectCores() - 1
}

registerDoParallel(cores = num_cores)

cat("Number of cores being used (1 free):", num_cores, "\n")
cat("SIMULATION: revision continuous \n")

M <- 1000
n_train <- 100
sigma <- sqrt(3)
N_post <- 5000
N_burn <- 2500

# RUN THE RESULTS HERE: three times p = 500, 1000, 50000
p_values <- c(500, 1000, 5000)  # adjust if needed

combined_results <- list()

run_simulation_tidy <- function(
    n_rep = 50, n = 100, p = 1000, sigma = 1, 
    N_post = 2000, N_burn = 2000, n_cores = 5
) {
  
  # foreach will now return a data frame for each 'i', and then rbind them together
  results_df <- foreach(
    i = 1:n_rep,
    .combine = 'rbind',
    .packages = c("ShrinkageTrees", "BART", "stochtree", "SoftBart")
  ) %dopar% {
    
    # Run the simulation for one iteration
    res_mat <- run_one_sim(seed = i, n = n, p = p, sigma = sigma, 
                           N_post = N_post, N_burn = N_burn)
    
    # Convert matrix to a data frame and add metadata
    tmp <- as.data.frame(res_mat)
    tmp$Method <- rownames(res_mat)
    tmp$Iter   <- i
    tmp$p      <- p
    
    tmp
  }
  
  return(results_df)
}

# --- Execution Logic ---
all_sim_data <- list()

for (p_val in p_values) {
  cat("Running simulation for p =", p_val, "\n")
  
  sim_df <- run_simulation_tidy(
    n_rep   = M,
    n       = n_train,
    p       = p_val,
    sigma   = sigma,
    N_post  = N_post,
    N_burn  = N_burn,
    n_cores = num_cores
  )
  
  all_sim_data[[as.character(p_val)]] <- sim_df
}

# Combine everything into one giant "Master Table"
final_flat_df <- do.call(rbind, all_sim_data)

# Define output file path 
# (! NAME MUST BE FILENAME_output.rds !)
output_file <- file.path(Sys.getenv("TMPDIR"), "revision_simple_output.rds")
cat("Saving all settings results to:", output_file, "\n")
saveRDS(final_flat_df, file = output_file)
cat("All results successfully saved in one file.\n")




library(dplyr)
final_flat_df <- readRDS("~/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/revision/revision_simple_output.rds")

summary_table <- final_flat_df %>%
  group_by(p, Method) %>%
  summarise(
    mean_RMSE     = mean(RMSE),
    mean_Coverage = mean(Coverage),
    mean_Length   = mean(Length),
    .groups = "drop"
  )

print(summary_table)


library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)

#--------------------------------------------------
# 1. Filter + Rename
#--------------------------------------------------

summary_table_filtered <- summary_table %>%
  filter(Method %in% c("CHF0.1", "BART", "DART", "SoftBART")) %>%
  mutate(
    Method = case_when(
      Method == "CHF0.1" ~ "Horseshoe Forest",
      TRUE ~ Method
    )
  )

df_long <- summary_table_filtered %>%
  rename(
    RMSE = mean_RMSE,
    Coverage = mean_Coverage,
    Length = mean_Length
  ) %>%
  pivot_longer(
    cols = c(RMSE, Coverage, Length),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = c("RMSE", "Coverage", "Length")),
    p = factor(p, levels = c(500, 1000, 5000)),
    Method = factor(
      Method,
      levels = c("Horseshoe Forest", "BART", "DART", "SoftBART")
    )
  )

#--------------------------------------------------
# 2. Plot
#--------------------------------------------------

p <- ggplot(df_long, aes(x = p, y = Value, fill = Method)) +

  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    color = "black",
    linewidth = 0.7
  ) +

  geom_hline(
    data = df_long %>% filter(Metric == "Coverage"),
    aes(yintercept = 0.95),
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +

  facet_grid2(
    . ~ Metric,
    scales = "free_y",
    independent = "y"
  ) +

  facetted_pos_scales(
    y = list(
      Metric == "RMSE" ~
        scale_y_continuous(
          limits = c(0.6, 1.2),
          breaks = c(0.7, 0.9, 1.1),
          expand = c(0, 0)
        ),

      Metric == "Coverage" ~
        scale_y_continuous(
          limits = c(0, 1.05),
          breaks = c(0.25, 0.5, 0.75, 1),
          expand = c(0, 0)
        ),

      Metric == "Length" ~
        scale_y_continuous(
          limits = c(0, 4),
          breaks = c(1, 2, 3),
          expand = c(0, 0)
        )
    )
  ) +

  labs(
    x = "Number of covariates",
    y = NULL,
    fill = "Method"
  ) +

  scale_fill_manual(values = c(
    "Horseshoe Forest" = "#042E78",  # deep blue
    "BART"             = "#FFF3E0",  # very pale cream
    "DART"             = "#FFCC80",  # soft apricot
    "SoftBART"         = "#FFB74D"   # orange
  )) +

  theme_bw(base_size = 21) +
  theme(
    strip.background = element_rect(fill = "white", color = NA),
    strip.text       = element_text(face = "bold"),
    axis.title.x     = element_text(face = "italic"),
    text             = element_text(family = "Times New Roman"),
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(p)

base_dir <- "~/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/revision"
ggsave(
  filename = file.path(base_dir, "revision_simple_bar_plot.png"),
  plot = p,
  width = 1.2 * 250,
  height = 0.9 * 150,
  units = "mm",
  dpi = 320
)









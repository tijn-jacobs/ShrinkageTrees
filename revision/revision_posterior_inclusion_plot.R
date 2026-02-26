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
    20 * (X[, 3] - 0.6)^2 +
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
  
  st_fit <- ShrinkageTrees::ShrinkageTrees(
    y = data$Y,
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
    f_hat  = as.numeric(st_fit$train_predictions),
    f_post = st_fit$train_predictions_sample,
    sigma_post = st_fit$sigma
  )
}


fit_bart <- function(
    data,
    N_post = 2500,
    N_burn = 2500
) {
  
  bart_fit <- ShrinkageTrees::SurvivalBART(
    X_train = data$X,
    time = data$Y,
    status = rep(1, length(data$Y)),
    timescale = "log",
    N_burn   = N_burn,
    N_post  = N_post,
    verbose = FALSE
  )
  
  list(
    f_hat  = as.numeric(bart_fit$train_predictions),
    f_post = bart_fit$train_predictions_sample,
    sigma_post = bart_fit$sigma
  )
}

fit_dart <- function(
    data,
    N_post = 2500,
    N_burn = 2500
) {
  
  dart_fit <- ShrinkageTrees::SurvivalDART(
    X_train = data$X,
    time = data$Y,
    status = rep(1, length(data$Y)),
    timescale = "log",
    N_burn   = N_burn,
    N_post  = N_post,
    verbose = FALSE
  )
  
  list(
    f_hat  = as.numeric(dart_fit$train_predictions),
    f_post = dart_fit$train_predictions_sample,
    sigma_post = dart_fit$sigma
  )
}


data <- data_gen(500, 10, sqrt(3))

dart_fit <- ShrinkageTrees::SurvivalDART(
  X_train = data$X,
  time = data$Y,
  status = rep(1, length(data$Y)),
  timescale = "log",
  N_burn   = 5000,
  N_post  = 5000,
  verbose = FALSE
)
dart_split_probs <- colMeans(dart_fit$split_probs)

bart_fit <- ShrinkageTrees::SurvivalBART(
  X_train = data$X,
  time = data$Y,
  status = rep(1, length(data$Y)),
  timescale = "log",
  N_burn   = 5000,
  N_post  = 5000,
  verbose = FALSE
)
M <- bart_fit$store_split_counts
M_norm <- M / rowSums(M)
bart_split_prob <- colMeans(M_norm)
bart_split_prob

m<-200
st_fit <- ShrinkageTrees::ShrinkageTrees(
  y = data$Y,
  X_train = data$X,
  outcome_type = "continuous",
  number_of_trees = m,
  prior_type = "horseshoe",
  local_hp  = k / sqrt(m),
  global_hp = k / sqrt(m),
  N_post = 5000,
  N_burn = 5000,
  store_posterior_sample = TRUE,
  verbose = FALSE
)
S <- st_fit$store_split_counts
S_norm <- S / rowSums(S)
st_split_prob <- colMeans(M_norm)
st_split_prob

par(mfrow=c(1,3))
plot(x=1:10, y=st_split_prob, type = "h", 
     ylim=c(0,0.6), 
     xlab = "Covariate index",
     ylab = "Posterior prob?",
     main = "Horseshoe Forest")
plot(x=1:10, y=bart_split_prob, type = "h",
     ylim=c(0,0.6), 
     xlab = "Covariate index",
     ylab = "Posterior prob?",
     main = "BART")
plot(x=1:10, y=dart_split_probs, type = "h", 
     ylim=c(0,0.6), 
     xlab = "Covariate index",
     ylab = "Posterior prob?",
     main = "DART")



















































set.seed(1)

# Number of repetitions
n_rep <- 3
p <- 10
n <- 500
N_post <- 5000
N_burn <- 2500

# Storage
all_st  <- vector("list", n_rep)
all_bart <- vector("list", n_rep)
all_dart <- vector("list", n_rep)

for (r in 1:n_rep) {
  
  cat("Running replication", r, "\n")
  
  data <- data_gen(n, p, sqrt(3))
  
  #------------------------------
  # 1) DART
  #------------------------------
  dart_fit <- ShrinkageTrees::SurvivalDART(
    X_train = data$X,
    time = data$Y,
    status = rep(1, length(data$Y)),
    timescale = "log",
    N_burn = N_burn,
    N_post = N_post,
    verbose = FALSE
  )
  
  all_dart[[r]] <- colMeans(dart_fit$split_probs)
  
  #------------------------------
  # 2) BART
  #------------------------------
  bart_fit <- ShrinkageTrees::SurvivalBART(
    X_train = data$X,
    time = data$Y,
    status = rep(1, length(data$Y)),
    timescale = "log",
    N_burn = N_burn,
    N_post = N_post,
    verbose = FALSE
  )
  
  M <- bart_fit$store_split_counts
  M_norm <- M / rowSums(M)
  all_bart[[r]] <- colMeans(M_norm)
  
  #------------------------------
  # 3) Horseshoe Forest
  #------------------------------
  k <- 0.1
  m <- 200
  
  st_fit <- ShrinkageTrees::ShrinkageTrees(
    y = data$Y,
    X_train = data$X,
    outcome_type = "continuous",
    number_of_trees = m,
    prior_type = "horseshoe",
    local_hp  = k / sqrt(m),
    global_hp = k / sqrt(m),
    N_post = N_post,
    N_burn = N_burn,
    verbose = FALSE
  )
  
  S <- st_fit$store_split_counts
  S_norm <- S / rowSums(S)
  all_st[[r]] <- colMeans(S_norm)
}

#--------------------------------------
# 3 x 3 Grid Plot
#--------------------------------------
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))

odd_breaks <- seq(1, p, by = 2)

for (r in 1:n_rep) {
  
  plot(1:p, all_st[[r]], type = "h",
       ylim = c(0, 0.6),
       xaxt = "n",
       xlab = "",
       ylab = "Posterior split prob",
       main = paste("Horseshoe Forest\nDataset", r))
  axis(1, at = odd_breaks)
  
  plot(1:p, all_bart[[r]], type = "h",
       ylim = c(0, 0.6),
       xaxt = "n",
       xlab = "",
       ylab = "Posterior split prob",
       main = paste("BART\nDataset", r))
  axis(1, at = odd_breaks)
  
  plot(1:p, all_dart[[r]], type = "h",
       ylim = c(0, 0.6),
       xaxt = "n",
       xlab = "",
       ylab = "Posterior split prob",
       main = paste("DART\nDataset", r))
  axis(1, at = odd_breaks)
}










library(dplyr)
library(tidyr)
library(ggplot2)

# Combine results into one long data frame
make_df <- function(probs_list, method_name) {
  bind_rows(lapply(seq_along(probs_list), function(r) {
    data.frame(
      Dataset = paste("Dataset", r),
      Variable = 1:length(probs_list[[r]]),
      SplitProb = probs_list[[r]],
      Method = method_name
    )
  }))
}

df_st   <- make_df(all_st,   "Horseshoe Forest")
df_bart <- make_df(all_bart, "BART")
df_dart <- make_df(all_dart, "DART")

df_all <- bind_rows(df_st, df_bart, df_dart)

# Color group
df_all <- df_all %>%
  mutate(Group = ifelse(Variable <= 5, "Signal", "Noise"))

# Plot
ggplot(df_all, aes(x = Variable, y = SplitProb)) +
  geom_segment(aes(xend = Variable,
                   y = 0,
                   yend = SplitProb,
                   color = Group),
               linewidth = 1.5) +
  facet_grid(Dataset ~ Method) +
  scale_color_manual(values = c(
    "Signal" = "chartreuse4",
    "Noise"  = "red3"
  )) +
  scale_x_continuous(breaks = seq(1, 10, by = 2)) +
  coord_cartesian(ylim = c(0, 0.65)) +
  labs(
    x = "Covariate index",
    y = "Posterior inclusion probability",
    color = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    text = element_text(family = "Times New Roman"),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

k<-12
ggsave(
  filename = "/Users/tijnjacobs/Documents/GitHub/ShrinkageTrees/revision/revision_posterior_inclusion_plot.png",
  width = k*185,
  height = k*100,
  units = "px",
  dpi = 320
)


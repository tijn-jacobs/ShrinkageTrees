### Simulation deeper Friedman ###

library(ShrinkageTrees)
library(AFTrees)
library(foreach)
library(doParallel)


FindKappa <- function(q, sigsq.hat, nu) {
  
  ### don't compute q directly if equals .9 or .99
  if(q < .905 & q > .895) {
    Q <- 6.1423
  }
  else if(q < .995 & q > .985) {
    Q <- 27.127
  }
  else if(q < .8 & q > .7) {
    Q <- 3.489662
  }
  else if(q < .55 & q > .45) {
    Q <- 2.287994
  }
  else if(q < .3 & q > .2) {
    Q <- 1.735869
  }
  else {
    Q <- qchiconv(q, nu=nu)
  }
  Kap <- sigsq.hat/Q
  return(Kap)
}

IndivAFT_new <- function(
    x.train, y.train, status, Trt, x.test=NULL,
    sigest=NA, sigdf=3, sigquant=.5,
    k=2.0,power=2.0, base=.95, nonparametric=TRUE,
    ntree=200,
    ndpost=1000, nskip=100,
    printevery=100, keepevery=1, keeptrainfits=TRUE,
    usequants=TRUE, numcut=100, printcutoffs=0,
    verbose=TRUE, scale="log")
{
  
  if(is.vector(x.train) | is.factor(x.train)) x.train = data.frame(x=x.train)
  
  if(!is.matrix(x.train)) {
    stop("x.train must be a matrix")
  }
  Trt.alt <- 1 - Trt ## counterfactual treatment
  if(is.null(x.test)) {
    x.test <- cbind(x.train, Trt.alt)
    x.train <- cbind(x.train, Trt)
    n.test <- 0
  } else {
    n.test <- nrow(x.test)
    t1 <- cbind(x.train, Trt.alt)
    t2 <- rbind(cbind(x.test, rep(1, n.test)), cbind(x.test, rep(0,n.test)))
    #print(dim(x.train))
    #print(length(Trt.alt))
    x.train <- cbind(x.train, Trt)
    x.test <- rbind(t1, t2)
  }
  
  #check input arguments:
  if((!is.matrix(x.train)) || (typeof(x.train)!="double")) stop("argument x.train must be a double matrix")
  if((!is.matrix(x.test)) || (typeof(x.test)!="double")) stop("argument x.test must be a double matrix")
  if((!is.vector(y.train)) || (typeof(y.train)!="double")) stop("argument y.train must be a double vector")
  
  if(nrow(x.train) != length(y.train)) stop("number of rows in x.train must equal length of y.train")
  if((nrow(x.test) >0) && (ncol(x.test)!=ncol(x.train))) stop("input x.test must have the same number of columns as x.train")
  if((!is.na(sigest)) && (typeof(sigest)!="double")) stop("input sigest must be double")
  if((!is.na(sigest)) && (sigest<0.0)) stop("input sigest must be positive")
  if((mode(sigdf)!="numeric") || (sigdf<0)) stop("input sigdf must be a positive number")
  if((mode(printevery)!="numeric") || (printevery<0)) stop("input printevery must be a positive number")
  if((mode(keepevery)!="numeric") || (keepevery<0)) stop("input keepevery must be a positive number")
  if((mode(sigquant)!="numeric") || (sigquant<0)) stop("input sigquant must be a positive number")
  if((mode(ntree)!="numeric") || (ntree<0)) stop("input ntree must be a positive number")
  if((mode(ndpost)!="numeric") || (ndpost<0)) stop("input ndpost must be a positive number")
  if((mode(nskip)!="numeric") || (nskip<0)) stop("input nskip must be a positive number")
  if((mode(k)!="numeric") || (k<0)) stop("input k must be a positive number")
  if(mode(numcut)!="numeric") stop("input numcut must be a numeric vector")
  if(length(numcut)==1) numcut = rep(numcut,ncol(x.train))
  if(length(numcut) != ncol(x.train)) stop("length of numcut must equal number of columns of x.train")
  numcut = as.integer(numcut)
  if(min(numcut)<1) stop("numcut must be >= 1")
  if(typeof(usequants)  != "logical") stop("input usequants must a logical variable")
  if(typeof(keeptrainfits)  != "logical") stop("input keeptrainfits must a logical variable")
  if(typeof(verbose)  != "logical") stop("input verbose must a logical variable")
  if(mode(printcutoffs)  != "numeric") stop("input printcutoffs must be numeric")
  printcutoffs = as.integer(printcutoffs)
  if(printcutoffs <0) stop("input printcutoffs must be >=0")
  if(power <= 0) stop("power must be positive")
  if(base <= 0) stop("base must be positive")
  
  null_lm <- survreg(Surv(y.train, status) ~ 1, dist="lognormal")
  null_sig <- null_lm$scale
  null_intercept <- null_lm$coefficients
  
  y_centered_log <- log(y.train) - null_intercept
  imr <- dnorm(log(y.train), mean=null_intercept, sd=null_sig)/pnorm(log(y.train), mean=null_intercept, sd=null_sig, lower.tail=FALSE)
  #y_unobs <- y.train*exp(.02*(1-status)*null_sig)   ### rough, estimate of unobserved survival times
  y_unobs <- exp( status*log(y.train) + (1 - status)*null_sig + (1 - status)*imr)
  
  #y.train.log = log(y.train)  ### transform to log-survival times
  rgy = range(log(y_unobs))
  zeta_tmp <- rgy[2] - rgy[1]
  zeta <- 4*null_sig
  #print(c(zeta,zeta_tmp))
  #aa = (rgy[1] + rgy[2])/(2*(rgy[2] - rgy[1]))
  
  
  # Need to include the survival package
  if (is.na(sigest)) {
    #templm = lm(y~x.train)
    #sigest = summary(templm)$sigma
    #tmplm <- survreg(Surv(y.train, status) ~ x.train, dist="lognormal")
    #sigest <- tmplm$scale
    sigest <- sd(log(y.train[status == 1]))
  } else {
    sigest = sigest #put input sigma estimate on transformed scale
  }
  
  ncskip = floor(nskip/keepevery)
  ncpost = floor(ndpost/keepevery)
  nctot = ncskip + ncpost
  totnd = keepevery*nctot
  
  nclust = 200
  kappa <- FindKappa(q=sigquant, sigsq.hat=sqrt(sigest), nu=sigdf)
  
  npind = ifelse(nonparametric, 1, 0)
  cres = .C('mbart',as.integer(nrow(x.train)), as.integer(ncol(x.train)), as.integer(nrow(x.test)),
            as.double(x.train), as.double(y_centered_log),
            as.double(x.test), as.integer(status),
            as.double(sigest),   as.integer(sigdf), as.double(sigquant),
            as.double(k),
            as.double(power), as.double(base),
            as.integer(ntree),      as.integer(totnd),
            as.integer(printevery), as.integer(keepevery),  as.integer(keeptrainfits),
            as.integer(numcut), as.integer(usequants), as.integer(printcutoffs),
            as.integer(verbose),
            sdraw=double(nctot),
            trdraw=double(nrow(x.train)*nctot),
            tedraw=double(nrow(x.test)*nctot),
            vcdraw=integer(ncol(x.train)*nctot), as.integer(nclust),
            mixdraw=double(nclust*nctot),
            locdraw=double(nclust*nctot),
            Mdraw=double(nctot),npind=as.integer(npind),as.double(kappa),
            as.double(zeta))
  # now read in the results...
  ## look at this: sigma is multiplied by (rgy[2] - rgy[1])
  sigma = cres$sdraw
  # sigma = cres$sdraw
  first.sigma = sigma[1:ncskip] # we often want the sigma draws
  sigma = sigma[ncskip+(1:ncpost)]
  
  # put sigest on the original y scale for output purposes
  sigest = sigest
  
  m.train = m.test = m.train.mean = m.test.mean = NULL
  varcount = NULL
  
  if (keeptrainfits) {
    m.train = matrix(cres$trdraw,nrow=nctot,byrow=T)[(ncskip+1):nctot,]
    m.train = m.train + null_intercept
    m.train.mean <- colMeans(m.train)
  }
  m.test = matrix(cres$tedraw,nrow=nctot,byrow=T)[(ncskip+1):nctot,]
  m.test = m.test + null_intercept
  
  mix.prop = matrix(cres$mixdraw,nrow=nctot,byrow=T)[(ncskip+1):nctot,]
  locations = matrix(cres$locdraw,nrow=nctot,byrow=T)[(ncskip+1):nctot,]
  mass = cres$Mdraw
  
  varcount = matrix(cres$vcdraw,nrow=nctot,byrow=T)[(ncskip+1):nctot,]
  
  npatients <- nrow(x.train)
  nsamps <- nrow(m.train)
  ThetaMat <- matrix(0, nrow=nsamps, ncol=npatients)
  Theta.test <- NULL
  
  if(scale=="log") {
    for(k in 1:npatients) {
      if(Trt[k]==1) {
        ThetaMat[,k] <- m.train[,k] - m.test[,k]
      }
      else {
        ThetaMat[,k] <- m.test[,k] - m.train[,k]
      }
    }
    if(n.test > 0) {
      Theta.test <- matrix(0, nrow=nsamps, ncol=n.test)
      for(h in 1:n.test) {
        Theta.test[,h] <- m.test[,npatients + h] - m.test[,npatients + h + n.test]
      }
    }
  }
  else if(scale=="time") {
    mgfs <- rowSums(exp(locations)*mix.prop)
    mgfs <- 1
    for(k in 1:npatients) {
      if(Trt[k]==1) {
        ThetaMat[,k] <- mgfs*(exp(m.train[,k] - m.test[,k]))
      }
      else {
        ThetaMat[,k] <- mgfs*(exp(m.test[,k] - m.train[,k]))
      }
    }
  }
  
  retval = list(
    call=match.call(),
    first.sigma=first.sigma,
    sigma=sigma,
    sigest=sigest,
    m.train=m.train,
    m.train.mean=m.train.mean,
    m.test=m.test,
    m.test.mean=m.test.mean,
    varcount=varcount,
    y = y.train,
    mix.prop=mix.prop,
    locations=locations,
    mass = mass,
    Theta = ThetaMat,
    Theta.test = Theta.test
  )
  class(retval) = 'indivaft'
  return(invisible(retval))
}


evaluate_IndivAFT <- function(data, use_uncensored = FALSE, ...) {
  if (!requireNamespace("AFTrees", quietly = TRUE)) {
    stop("Package 'AFTrees' is required.")
  }
  
  # Determine which outcome to use
  if (use_uncensored) {
    y_train <- exp(data$uncensored_event_times)
    status <- rep(1, length(data$status))
  } else {
    y_train <- exp(data$follow_up)
    status <- data$status
  }
  
  # Suppress ALL output from IndivAFT
  # Call IndivAFT and suppress all output correctly
  fit <- IndivAFT_new(
    x.train = data$X_train,
    y.train = y_train,
    status = status,
    Trt = data$treatment,
    x.test = data$X_train,
    sigdf = 3,
    sigquant = 0.5,
    k = 2.0,
    power = 2.0,
    base = 0.95,
    printevery = 0,
    keepevery = 1,
    keeptrainfits = TRUE,
    nonparametric = FALSE,
    scale = "time",
    verbose = FALSE,
    ...
  )
  
  # Get log-transformed CATE posterior samples
  fit$Theta <- log(fit$Theta)
  
  # ATE
  ate_samples <- rowMeans(fit$Theta)
  ate_estimate <- mean(ate_samples)
  ate_ci <- quantile(ate_samples, probs = c(0.025, 0.975))
  ate_coverage <- unname(ate_ci[1] <= data$true_ate & ate_ci[2] >= data$true_ate)
  ate_ci_length <- unname(diff(ate_ci))
  
  # CATE
  cate_mean <- colMeans(fit$Theta)
  cate_rpehe <- sqrt(mean((cate_mean - data$true_cate)^2))
  
  cate_ci <- apply(fit$Theta, 2, quantile, probs = c(0.025, 0.975))
  cate_coverage <- mean(cate_ci[1, ] <= data$true_cate & cate_ci[2, ] >= data$true_cate)
  cate_ci_length <- mean(cate_ci[2, ] - cate_ci[1, ])
  
  # RMSE for prediction of survival times
  m_train_log <- fit$m.train.mean
  rmse <- sqrt(mean((m_train_log - data$true_event_times)^2))
  
  # Compute the C-index
  if (requireNamespace("survival", quietly = TRUE)) {
    c_index <- survival::concordance(
      survival::Surv(exp(data$follow_up), data$status) ~ exp(m_train_log)
    )$concordance
  } else {
    c_index <- NA
  }
  
  return(list(
    postmean_sigma = mean(fit$sigma),
    ate = ate_estimate,
    ate_ci = ate_ci,
    ate_coverage = ate_coverage,
    ate_ci_length = ate_ci_length,
    cate_rpehe = cate_rpehe,
    cate_coverage = cate_coverage,
    cate_ci_length = cate_ci_length,
    rmse = rmse,
    c_index = c_index
  ))
}


evaluate_CHF <- function(data, use_uncensored = FALSE, ...) {
  if (!requireNamespace("ShrinkageTrees", quietly = TRUE)) {
    stop("Package 'ShrinkageTrees' is required.")
  }
  
  # Determine outcome and censoring indicator
  if (use_uncensored) {
    y <- data$uncensored_event_times
    status <- NULL
    scale <- "log"
  } else {
    y <- data$follow_up
    status <- data$status
    scale <- "log"
  }
  
  X1 <- cbind(data$propensity, data$X_train)

  # Fit the CHF model
  fit <- ShrinkageTrees::CausalShrinkageForest(
    y = y,
    status = status,
    outcome_type = "right-censored",
    X_train_control = X1,
    X_train_treat = data$X_train,
    prior_type_control = "horseshoe",
    prior_type_treat = "horseshoe",
    treatment_indicator_train = data$treatment,
    store_posterior_sample = TRUE,
    timescale = scale,
    verbose = FALSE,
    ...
  )
  
  # Posterior samples of τ(x)
  cate_samples <- fit$train_predictions_sample_treat
  
  # Pointwise quantities
  cate_mean <- fit$train_predictions_treat
  cate_rpehe <- sqrt(mean((cate_mean - data$true_cate)^2))
  
  # Per-observation posterior intervals for τ(x)
  cate_ci <- apply(cate_samples, 2, quantile, probs = c(0.025, 0.975))
  cate_coverage <- mean(cate_ci[1,] <= data$true_cate & cate_ci[2,] >= data$true_cate)
  cate_ci_length <- mean(cate_ci[2,] - cate_ci[1,])
  
  # ATE estimate and uncertainty
  ate_samples <- rowMeans(cate_samples)
  ate_estimate <- mean(ate_samples)
  ate_ci <- quantile(ate_samples, probs = c(0.025, 0.975))
  ate_coverage <- unname(ate_ci[1] <= data$true_ate & ate_ci[2] >= data$true_ate)
  ate_ci_length <- unname(diff(ate_ci))
  
  # RMSE for total outcome prediction
  rmse <- sqrt(mean((fit$train_predictions - data$true_event_times)^2))
  
  # Compute the C-index
  if (requireNamespace("survival", quietly = TRUE)) {
    c_index <- survival::concordance(
      survival::Surv(exp(data$follow_up), data$status) ~ exp(fit$train_predictions)
    )$concordance
  } else {
    c_index <- NA
  }
  
  return(list(
    postmean_sigma = mean(fit$sigma),
    ate = ate_estimate,
    ate_ci = ate_ci,
    ate_coverage = ate_coverage,
    ate_ci_length = ate_ci_length,
    cate_rpehe = cate_rpehe,
    cate_coverage = cate_coverage,
    cate_ci_length = cate_ci_length,
    rmse = rmse,
    c_index = c_index
  ))
}


data_gen <- function(n_train, p_feat, sigma, s_prog, cens_scale = 1) {
  X_train <- matrix(runif(n_train * p_feat), n_train, p_feat)
  
  propensity <- pnorm(-0.5 + 0.4*X_train[,1])
  treatment <- rbinom(n_train, 1, propensity)
  
  beta_prog <- rnorm(p_feat, 0, 1) * rbinom(p_feat, 1, s_prog)
  mu <- beta_prog %*% t(X_train)
  
  tau <- 10 * sin(pi * X_train[, 1] * X_train[, 2]) +
    20 * (X_train[, 3] - 0.5)^2 +
    10 * X_train[, 4] +
    5  * X_train[, 5]
  
  true_ate <- 10 * 0.5246745 + 20 * 1/12 + 10 * 1/2 + 5 * 1/2
  
  true_event_times <- mu + (treatment - 0.5) * tau
  uncensored_event_times <- true_event_times + rnorm(n_train, 0, sigma)
  sd_un <- sd(uncensored_event_times)
  uncensored_event_times <- uncensored_event_times / sd_un
  true_event_times <- true_event_times / sd_un
  
  C <- log(rexp(n_train, cens_scale)) + min(uncensored_event_times)
  follow_up <- pmin(uncensored_event_times, C)
  status <- as.numeric(uncensored_event_times <= C)
  
  return(list(
    X_train = X_train,
    treatment = treatment,
    propensity = propensity,
    follow_up = as.numeric(follow_up),
    status = status,
    true_event_times = as.numeric(true_event_times),
    uncensored_event_times = as.numeric(uncensored_event_times),
    true_cate = as.vector(tau) / sd_un,
    true_ate = true_ate / sd_un,
    sample_ate = mean(tau) / sd_un,
    obs_sigma = sigma / sd_un
  ))
}


run_single_simulation <- function(n_train,
                                  p_feat,
                                  sigma,
                                  s_prog,
                                  cens_scale,
                                  k,
                                  omega,
                                  N_post = 2000,
                                  N_burn = 2000,
                                  seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Generate data
  dt <- data_gen(n_train, p_feat, sigma, s_prog, cens_scale)
  
  # Evaluate IndivAFT
  results_indiv_aft <- evaluate_IndivAFT(
    dt,
    nskip = N_burn,
    ndpost = N_post,
    ntree = 200
  )
  
  # Evaluate CHF with k = 0.1
  results_chf0.1 <- evaluate_CHF(
    dt,
    number_of_trees_control = 200,
    number_of_trees_treat = 200,
    local_hp_control  = 0.1 / sqrt(200),
    global_hp_control = 0.1 / sqrt(200),
    local_hp_treat    = 0.1 / sqrt(200),
    global_hp_treat   = 0.1 / sqrt(200),
    N_post = N_post,
    N_burn = N_burn
  )
  
  # Return combined results
  return(data.frame(
    Method = c("IndivAFT", "CHF (k=0.1)"),
    sigma = c(results_indiv_aft$postmean_sigma,
              results_chf0.1$postmean_sigma),
    ATE = c(results_indiv_aft$ate,
            results_chf0.1$ate),
    ATE_bias = c(results_indiv_aft$ate - dt$true_ate,
                 results_chf0.1$ate - dt$true_ate),
    ATE_coverage = c(results_indiv_aft$ate_coverage,
                     results_chf0.1$ate_coverage),
    ATE_CI_Length = c(results_indiv_aft$ate_ci_length,
                      results_chf0.1$ate_ci_length),
    CATE_RPEHE = c(results_indiv_aft$cate_rpehe,
                   results_chf0.1$cate_rpehe),
    CATE_coverage = c(results_indiv_aft$cate_coverage,
                      results_chf0.1$cate_coverage),
    CATE_CI_Length = c(results_indiv_aft$cate_ci_length,
                       results_chf0.1$cate_ci_length),
    RMSE = c(results_indiv_aft$rmse,
             results_chf0.1$rmse),
    C_Index = c(results_indiv_aft$c_index,
                results_chf0.1$c_index)
  ))
}






run_parallel_simulations <- function(M,
                                     n_train,
                                     p_feat,
                                     sigma,
                                     s_prog,
                                     cens_scale,
                                     k,
                                     omega,
                                     N_post = 2000,
                                     N_burn = 2000) {
  
  results <- foreach(
    i = 1:M,
    .combine = rbind,
    .packages = c("ShrinkageTrees", "AFTrees", "survival"),
    .export = c("FindKappa", "IndivAFT_new", "evaluate_IndivAFT", "evaluate_CHF", 
                "data_gen", "run_single_simulation")  # list all needed local functions
  ) %dopar% {
    run_single_simulation(
      n_train = n_train,
      p_feat = p_feat,
      sigma = sigma,
      s_prog = s_prog,
      cens_scale = cens_scale,
      k = k,
      omega = omega,
      N_post = N_post,
      N_burn = N_burn,
      seed = i
    )
  }
  
  return(results)
}

run_simulations_over_p <- function(
    p_values,
    cens_scales,
    M,
    n_train,
    sigma,
    s_prog,
    k,
    omega,
    N_post = 2000,
    N_burn = 2000
) {
  if (length(p_values) != length(cens_scales)) {
    stop("p_values and cens_scales must be the same length.")
  }
  
  all_results <- Map(function(p_val, cs_val) {
    cat("Running for p =", p_val, "\n")
    res <- run_parallel_simulations(
      M = M,
      n_train = n_train,
      p_feat = p_val,
      sigma = sigma,
      s_prog = s_prog,
      cens_scale = cs_val,
      k = k,
      omega = omega,
      N_post = N_post,
      N_burn = N_burn
    )
    res$p_feat <- p_val  # retain p
    return(res)          # cens_scale not added
  }, p_values, cens_scales)
  
  do.call(rbind, all_results)
}

find_cens_scale <- function(p,
                            n_train = 100,
                            sigma = sqrt(3),
                            s_prog = 0.1,
                            target_cens,
                            M = 10000,
                            verbose = FALSE) {
  
  # Inner function: return absolute difference between actual and target censoring
  estimate_cens_diff <- function(cens_scale) {
    rates <- replicate(M, {
      dt <- data_gen(n_train, p, sigma, s_prog, cens_scale)
      mean(1 - dt$status)
    })
    diff <- abs(mean(rates) - target_cens)
    if (verbose) cat("cens_scale =", round(cens_scale, 4), 
                     "| estimated cens =", round(mean(rates), 3), "\n")
    return(diff)
  }
  
  # Optimize over a reasonable range
  opt_result <- optimize(estimate_cens_diff, interval = c(0.001, 2), tol = 1e-3)
  
  return(opt_result$minimum)
}


summarize_simulation_results_base <- function(results_df, digits = 3) {
  unique_methods <- unique(results_df$Method)
  unique_p <- sort(unique(results_df$p_feat))
  
  all_results <- list()
  
  for (p in unique_p) {
    for (method in unique_methods) {
      subset_df <- results_df[results_df$p_feat == p & results_df$Method == method, ]
      
      # Compute average metrics
      avg_metrics <- colMeans(subset_df[, c(
        "ATE_coverage", "ATE_CI_Length",
        "CATE_RPEHE", "CATE_coverage", "CATE_CI_Length",
        "RMSE", "C_Index"
      )], na.rm = TRUE)
      
      # Compute ATE RMSE separately
      ate_rmse <- sqrt(mean(subset_df$ATE_bias^2, na.rm = TRUE))
      
      all_results[[length(all_results) + 1]] <- c(
        p_feat = p,
        Method = method,
        ATE_RMSE = round(ate_rmse, digits),
        round(avg_metrics, digits)
      )
    }
  }
  
  summary_df <- do.call(rbind, all_results)
  rownames(summary_df) <- NULL
  return(as.data.frame(summary_df))
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
cat("SIMULATION: DEEPER FRIEDMAN2 high\n")

M <- 1000
n_train <- 100
sigma <- sqrt(3)
s_prog <- 0.1
N_post <- 5000
N_burn <- 2500


p_vals <- round(10^seq(log10(10), log10(1000), length.out = 20))
# cens_scales_high <- sapply(p_vals, function(p) find_cens_scale(p, target_cens = 0.60, M = 500))
cens_scales_high <- c(0.1948321, 0.1812497, 0.1951909, 0.1957655, 0.1963035, 0.1942478, 0.1933557, 0.1910802, 0.1926626, 0.1915547, 0.1884151, 0.1857072, 0.1830651, 0.1794495, 0.1812497, 0.1746715, 0.1699665, 0.167341, 0.1583239, 0.1558305)

# Do the HD scenario now
cens_scales_high <- cens_scales_high[p_vals <= 100]
p_vals <- p_vals[p_vals <= 100]


cat("\nStarting with high censoring scenario\n\n")
results_high <- run_simulations_over_p(
  p_values = p_vals,
  cens_scales = cens_scales_high,
  M = M,
  n_train = n_train,
  sigma = sigma,
  s_prog = s_prog,
  k = 0.1, # Redundant
  omega = 1/2,
  N_post = N_post,
  N_burn = N_burn
)


combined_results <- list(results_high = results_high)

# Define output file path 
# (! NAME MUST BE FILENAME_output.rds !)
output_file <- file.path(Sys.getenv("TMPDIR"), "sim_deeper_friedman2_high_LD_output.rds")

# Print message
cat("Saving all settings results to:", output_file, "\n")

# Save the combined list
saveRDS(combined_results, file = output_file)

# Confirm successful save
cat("All results successfully saved in one file.\n")


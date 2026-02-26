#' Compute mean estimate for censored data
#'
#' Estimates the mean and standard deviation for right-censored survival data.
#' Uses the `afthd` package if available (placeholder), else `survival`, and 
#' otherwise falls back to the naive mean among observed events.
#'
#' @param y Numeric vector of (log-transformed) survival times.
#' @param status Numeric vector; event indicator (1 = event, 0 = censored).
#'
#' @return A list with elements:
#'   \item{mu}{Estimated mean of survival times.}
#'   \item{sd}{Estimated standard deviation of survival times.}
#'   \item{min}{Estimated minimum of survival times.}
#'   \item{max}{Estimated maximum of survival times.}
#' 
#' @importFrom stats dnorm pnorm
censored_info <- function(y, status) {
  
  if (requireNamespace("survival", quietly = TRUE)) {
    
    # Fit Gaussian AFT with censoring
    fit <- survival::survreg(survival::Surv(y, status) ~ 1, dist = "gaussian")
    
    mu <- as.numeric(fit$coefficients[1])
    sd <- as.numeric(fit$scale)
    
    # ---- IMPUTE censored values ----
    cens_idx <- which(status == 0)
    imputed_y <- y
    
    if (length(cens_idx) > 0) {
      a <- (y[cens_idx] - mu) / sd
      lambda <- dnorm(a) / (1 - pnorm(a))     # Mills ratio
      imputed_y[cens_idx] <- mu + sd * lambda # Replace censored with E[Y | Y > c]
    }
    
    # min/max after imputation
    min_val <- min(imputed_y, na.rm = TRUE)
    max_val <- max(imputed_y, na.rm = TRUE)
    
  } else {
    
    # Fallback: use only uncensored observations
    min_val <- min(y, na.rm = TRUE)
    max_val <- max(y, na.rm = TRUE)
    
    mu <- mean(y[status == 1], na.rm = TRUE)
    sd <- sd(y[status == 1], na.rm = TRUE)
  }
  
  return(list(
    mu = mu,
    sd = sd,
    min = min_val,
    max = max_val
  ))
}

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
#'
censored_info <- function(y, status) {
  
  # Currently, version of afthd does not work
  # if (requireNamespace("afthd", quietly = TRUE)) {
  #   # Placeholder: afthd implementation can be added here later
  #   warning("afthd method not yet implemented; falling back to survreg or naive estimate.")
    
  # } else 
  
  if (requireNamespace("survival", quietly = TRUE)) {
    
    fit <- survival::survreg(survival::Surv(y, status) ~ 1, dist = "gaussian")
    mu <- as.numeric(fit$coefficients[1])
    sd <- as.numeric(fit$scale)
    
  } else {
    # Fallback: use mean and SD of uncensored observations
    mu <- mean(y[status == 1], na.rm = TRUE)
    sd <- sd(y[status == 1], na.rm = TRUE)
  }
  
  return(list(
    mu = mu,
    sd = sd
  ))
}

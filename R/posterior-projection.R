# Lower-dimensional posterior projection (Woody, Carvalho & Murray, 2021).
# Projects every posterior draw of a fitted function onto a simpler model;
# the user chooses the covariates, the method does no selection of its own.

#' Is this fit a survival model?
#'
#' `timescale` governs the stored scale only for survival outcomes. It keeps
#' its default (`"time"`) for continuous and binary fits, whose draws are never
#' exponentiated, so it must not be read on its own.
#'
#' @noRd
.pp_is_survival <- function(object)
  isTRUE(object$outcome_type %in% c("right-censored", "interval-censored"))

#' Extract posterior draws on the additive scale
#'
#' Returns the posterior draws of the fitted function (S x n) together with the
#' evaluation design (n x p). The stored draws are exponentiated when
#' `timescale = "time"`, and binary fits store probabilities, so neither can be
#' projected as stored; this function normalises both to the additive scale and
#' reports which scale it ended up on.
#'
#' @param object A fitted `ShrinkageTrees` or `CausalShrinkageForest` object.
#' @param target For causal fits, `"treatment"` (project tau) or
#'   `"prognostic"` (project mu). Ignored otherwise.
#' @return A list with `draws`, `X`, `scale` and `target`.
#' @noRd
.pp_extract <- function(object, target = c("treatment", "prognostic")) {

  target <- match.arg(target)

  if (inherits(object, "CausalShrinkageForest")) {
    if (target == "treatment") {
      draws <- object$train_predictions_sample_treat
      X     <- object$data$X_train_treat
    } else {
      draws <- object$train_predictions_sample_control
      X     <- object$data$X_train_control
    }
    if (is.null(draws))
      stop("No posterior sample stored. Refit with store_posterior_sample = TRUE.")
    ## timescale only governs the stored scale for survival outcomes; it keeps
    ## its default value for continuous fits, whose draws are never exponentiated.
    if (.pp_is_survival(object) && identical(object$timescale, "time")) {
      draws     <- log(draws)
      scale_lab <- "log-time (inverted from stored acceleration factors)"
    } else if (.pp_is_survival(object)) {
      scale_lab <- "log-time"
    } else {
      scale_lab <- "response"
    }

  } else if (inherits(object, "ShrinkageTrees")) {
    draws <- object$train_predictions_sample
    X     <- object$data$X_train
    if (is.null(draws))
      stop("No posterior sample stored. Refit with store_posterior_sample = TRUE.")
    if (identical(object$outcome_type, "binary")) {
      scale_lab <- "probit (latent)"
    } else if (!.pp_is_survival(object)) {
      scale_lab <- "response"
    } else if (identical(object$timescale, "time")) {
      draws     <- log(draws)
      scale_lab <- "log-time (inverted from stored time scale)"
    } else {
      scale_lab <- "log-time"
    }

  } else {
    stop(".pp_extract() expects a ShrinkageTrees or CausalShrinkageForest ",
         "object.")
  }

  X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  if (ncol(draws) != nrow(X))
    stop("draws has ", ncol(draws), " columns but the design has ", nrow(X),
         " rows -- check orientation.")

  list(draws = draws, X = X, scale = scale_lab, target = target)
}


.pp_wls <- function(B, Ft, w, ridge = 1e-8) {
  sw  <- sqrt(w)
  Bw  <- B * sw
  BtB <- crossprod(Bw)
  d   <- mean(diag(BtB))
  if (!is.finite(d) || d <= 0) d <- 1
  A    <- BtB + diag(ridge * d, ncol(Bw))
  coef <- solve(A, crossprod(Bw, Ft * sw))
  list(coef = coef, fitted = B %*% coef)
}

#' Choose lambda by cross-validation on the posterior mean
#'
#' `lambda` indexes the summary family, so it is fixed once, for all draws,
#' before any projection happens. `cv.glmnet` is run on \eqn{\bar f}: holding out
#' observations asks how well a sparse linear summary generalises across the
#' covariate space, which is a real question even though a posterior draw
#' carries no observation noise -- with `p >= n` the lasso can fit the training
#' points using covariates that predict nothing at held-out `x`.
#'
#' @param lambda `NULL` or `"1se"` for `lambda.1se`, `"min"` for `lambda.min`,
#'   or a number, which skips cross-validation.
#' @noRd
.pp_cv_lambda <- function(fbar, X, w, alpha, lambda) {
  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("penalised projections need the glmnet package.", call. = FALSE)
  if (is.numeric(lambda))
    return(list(lambda = lambda, rule = "user", cv = NULL))

  rule <- if (is.null(lambda)) "1se" else match.arg(lambda, c("1se", "min"))
  cv   <- glmnet::cv.glmnet(X, fbar, alpha = alpha, weights = w,
                            standardize = TRUE)
  list(lambda = if (rule == "min") cv$lambda.min else cv$lambda.1se,
       rule = rule, cv = cv)
}

#' Solve the penalised projection separately for every draw
#'
#' For each posterior draw,
#' \eqn{\gamma^{(s)} = \arg\min \|f^{(s)} - X\gamma\|^2 + \lambda\|\gamma\|_1}.
#' The lasso is a nonlinear operator, so the projection of the posterior mean is
#' *not* the mean of the projections; solving per draw is what makes the result
#' a posterior over sparse projections, with genuine selection uncertainty and
#' coefficients that are exactly zero in some draws and not in others.
#'
#' Costs one glmnet solve per draw. A short warm-start sequence ending at
#' `lambda` is used rather than a single value, which glmnet handles far better.
#'
#' @noRd
.pp_penalised_draws <- function(Ft, X, w, alpha, lambda, verbose = TRUE) {
  S    <- ncol(Ft)
  lseq <- lambda * c(8, 4, 2, 1)
  B    <- matrix(0, ncol(X) + 1L, S)

  if (verbose && S > 2000)
    message("Solving ", S, " penalised projections, one per draw ...")

  ## Take the last column of the path directly. coef(g, s = lambda)
  ## interpolates between neighbouring lambdas, which can blend a zero with a
  ## nonzero and destroy the exact sparsity we are here to report.
  k <- length(lseq)
  for (i in seq_len(S)) {
    g <- glmnet::glmnet(X, Ft[, i], alpha = alpha, lambda = lseq,
                        weights = w, standardize = TRUE)
    B[, i] <- c(g$a0[k], as.matrix(g$beta)[, k])
  }
  rownames(B) <- c("(Intercept)", colnames(X))

  ## glmnet can fail numerically on individual draws; a single non-finite
  ## column would otherwise poison every downstream summary (residuals,
  ## rho^2, sigma) via rowMeans. Drop such draws loudly, never silently.
  ok <- colSums(!is.finite(B)) == 0L
  if (!all(ok)) {
    warning(sum(!ok), " of ", S, " penalised projections returned ",
            "non-finite coefficients and were dropped from the summary.",
            call. = FALSE)
    B <- B[, ok, drop = FALSE]
  }

  fitted  <- X %*% B[-1L, , drop = FALSE] + rep(B[1L, ], each = nrow(X))
  nonzero <- rowMeans(B[-1L, , drop = FALSE] != 0)
  keep    <- c(TRUE, nonzero > 0)      # intercept plus anything ever selected

  list(coef    = B[keep, , drop = FALSE],
       fitted  = fitted,
       nonzero = c(1, nonzero)[keep],
       edf     = mean(colSums(B[-1L, , drop = FALSE] != 0)) + 1,
       kept    = which(ok))
}

.pp_rho2 <- function(Ft, fitted, w) {
  n   <- nrow(Ft)
  wm  <- colSums(w * Ft) / sum(w)
  tss <- colSums(w * (Ft - rep(wm, each = n))^2)
  rss <- colSums(w * (Ft - fitted)^2)
  ifelse(tss > 0, 1 - rss / tss, NA_real_)
}


.pp_basis_linear <- function(X) {
  B <- cbind(1, X)
  colnames(B) <- c("(Intercept)", colnames(X))
  B
}

.pp_basis_additive <- function(X, df_spline = 4) {
  parts <- list(); objs <- list(); idx <- list()
  pos <- 1L
  for (j in seq_len(ncol(X))) {
    xj <- X[, j]; nm <- colnames(X)[j]; u <- length(unique(xj))
    if (u <= 2) {
      m <- matrix(xj, ncol = 1, dimnames = list(NULL, nm)); ob <- NULL
    } else {
      ob <- splines::ns(xj, df = min(df_spline, u - 1L))
      m  <- as.matrix(ob)
      colnames(m) <- paste0(nm, ".s", seq_len(ncol(m)))
    }
    parts[[length(parts) + 1L]] <- m
    objs[[nm]] <- ob
    idx[[nm]]  <- pos + seq_len(ncol(m))
    pos <- pos + ncol(m)
  }
  B <- cbind(matrix(1, nrow(X), 1, dimnames = list(NULL, "(Intercept)")),
             do.call(cbind, parts))
  list(B = B, splines = objs, index = idx)
}

.pp_tree_structure <- function(f, X, w, max_leaves) {
  if (!requireNamespace("rpart", quietly = TRUE))
    stop("family = 'tree' needs the rpart package.")
  dat  <- data.frame(.f = f, as.data.frame(X), check.names = TRUE)
  ctrl <- rpart::rpart.control(cp = 0, xval = 0, maxcompete = 0,
                               maxsurrogate = 0,
                               minbucket = max(5L, round(0.02 * nrow(X))))
  tr  <- rpart::rpart(.f ~ ., data = dat, weights = w, control = ctrl)
  cpt <- as.data.frame(tr$cptable)

  ok <- cpt$nsplit[cpt$nsplit >= 1 & (cpt$nsplit + 1L) <= max_leaves]
  if (!length(ok))
    stop("rpart produced no split with at most ", max_leaves, " leaves.")
  cp_k <- cpt$CP[cpt$nsplit == max(ok)][1]
  pr   <- rpart::prune(tr, cp = cp_k * 1.0000001)

  leaf <- factor(pr$where)
  if (nlevels(leaf) < 2) stop("The pruned tree has a single leaf.")
  B <- stats::model.matrix(~ leaf - 1)
  colnames(B) <- paste0("leaf", seq_len(ncol(B)))

  list(B = B, size = nlevels(leaf), tree = pr,
       covariates = unique(as.character(pr$frame$var[pr$frame$var != "<leaf>"])),
       leaf_rows = sort(unique(pr$where)))
}

#' Check that the projection is defined
#'
#' The unpenalised projection is an ordinary least squares fit, so it needs
#' fewer basis columns than observations. With `p >= n` -- the regime this
#' package targets -- a linear projection onto every covariate is undefined,
#' and one onto anything approaching `n` covariates interpolates the draws and
#' returns a summary \eqn{R^2} near 1 whatever the model is doing.
#'
#' @noRd
.pp_check_dim <- function(q, n, family, p_used) {
  if (q >= n)
    stop("The projection needs ", q, " basis columns but there are only ", n,
         " observations, so it is undefined.\n",
         "  Either pass `covariates = ` a subset of at most ",
         if (family == "additive") "roughly " else "", n - 2L,
         " basis columns' worth,\n",
         "  or use penalty = \"ridge\" to project onto all ", p_used,
         " covariates with regularisation.", call. = FALSE)
  if (q > n / 10)
    warning("The summary uses ", q, " basis columns for ", n,
            " observations. The summary R^2 is inflated by the size of ",
            "the summary itself, not only by fit.", call. = FALSE)
  invisible(TRUE)
}

#' Posterior of each additive partial effect, 10th to 90th percentile
#'
#' The raw spline basis coefficients (`age.s1`, `age.s2`, ...) mean nothing on
#' their own. This reports, per covariate, how far its partial effect moves
#' between the 10th and 90th percentile of that covariate -- a nonlinear
#' analogue of a slope, on the same scale as the response.
#'
#' @noRd
.pp_additive_effects <- function(x, a) {
  ad <- x$structure
  qs <- vapply(names(ad$index), function(nm)
    stats::quantile(x$X[, nm], c(0.1, 0.9), names = FALSE), numeric(2))

  d <- vapply(seq_along(ad$index), function(i) {
    nm <- names(ad$index)[i]
    ob <- ad$splines[[nm]]
    B  <- if (is.null(ob)) matrix(qs[, i], ncol = 1) else
      stats::predict(ob, qs[, i])
    cf <- x$coefficients[, ad$index[[nm]], drop = FALSE]     # S x q
    as.numeric(cf %*% (B[2, ] - B[1, ]))
  }, numeric(nrow(x$coefficients)))
  colnames(d) <- names(ad$index)

  data.frame(
    variable = colnames(d),
    mean     = colMeans(d),
    sd       = apply(d, 2, sd),
    lower    = apply(d, 2, quantile, a),
    upper    = apply(d, 2, quantile, 1 - a),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

.pp_coef_summary <- function(coefs, a) {
  data.frame(
    variable = colnames(coefs),
    mean     = colMeans(coefs),
    sd       = apply(coefs, 2, sd),
    lower    = apply(coefs, 2, quantile, a),
    upper    = apply(coefs, 2, quantile, 1 - a),
    row.names = NULL, stringsAsFactors = FALSE
  )
}


#' Lower-dimensional posterior projection
#'
#' Projects every posterior draw of the fitted function onto a simpler model,
#' following Woody, Carvalho and Murray (2021), and returns the resulting
#' posterior distribution over projections.
#'
#' The projection is taken over the design the model was fitted to; there is no
#' `newdata` argument. Nor is there any covariate selection: by default the
#' projection uses every covariate, and a subset can be supplied through
#' `covariates`. Choosing that subset is the user's job.
#'
#' @references Woody, S., Carvalho, C. M., & Murray, J. S. (2021). Model
#'   Interpretation Through Lower-Dimensional Posterior Summarization.
#'   \emph{Journal of Computational and Graphical Statistics}, 30(1), 144-161.
#' @seealso [print.PosteriorProjection()], [plot.PosteriorProjection()]
#' @examples
#' n <- 80; p <- 30
#' X <- matrix(rnorm(n * p), n, p,
#'             dimnames = list(NULL, paste0("x", 1:p)))
#' y <- 2 * X[, 1] - X[, 2] + rnorm(n)
#'
#' fit <- HorseTrees(y = y, X_train = X, outcome_type = "continuous",
#'                   number_of_trees = 5, N_post = 50, N_burn = 25,
#'                   store_posterior_sample = TRUE, verbose = FALSE)
#'
#' # Linear projection onto two covariates
#' pp <- posterior_projection(fit, family = "linear",
#'                            covariates = c("x1", "x2"))
#' pp
#'
#' # Additive projection with per-covariate effect curves
#' pp_add <- posterior_projection(fit, family = "additive",
#'                                covariates = c("x1", "x2"))
#'
#' \donttest{
#' # Sparse projection; lambda cross-validated on the posterior mean
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   pp_las <- posterior_projection(fit, family = "linear",
#'                                  penalty = "lasso")
#' }
#'
#' # Subgroup summary
#' if (requireNamespace("rpart", quietly = TRUE)) {
#'   pp_tree <- posterior_projection(fit, family = "tree", max_leaves = 4)
#' }
#' }
#' @rdname posterior_projection
#' @export
posterior_projection <- function(object, ...) UseMethod("posterior_projection")

#' @rdname posterior_projection
#' @export
posterior_projection.ShrinkageTrees <- function(object, ...) {
  ex <- .pp_extract(object)
  posterior_projection.default(ex$draws, X = ex$X, scale_label = ex$scale,
                               target = "f", ...)
}

#' @param target Which function to project: `"treatment"` for the CATE surface
#'   \eqn{\tau(x)}, or `"prognostic"` for \eqn{\mu(x)}.
#' @rdname posterior_projection
#' @export
posterior_projection.CausalShrinkageForest <- function(
    object, target = c("treatment", "prognostic"), ...) {
  target <- match.arg(target)
  ex <- .pp_extract(object, target = target)
  posterior_projection.default(ex$draws, X = ex$X, scale_label = ex$scale,
                               target = ex$target, ...)
}

#' @param object A fitted `ShrinkageTrees` or `CausalShrinkageForest` object,
#'   or (default method) an `S` x `n` matrix of posterior draws on the
#'   additive scale.
#' @param X (default method) The `n` x `p` design the supplied draws
#'   correspond to; column names label the projection. Not a `newdata`
#'   argument -- it must come from somewhere when there is no fitted object
#'   to read it from.
#' @param family `"linear"` projects onto a linear model, `"additive"` onto a
#'   natural-spline additive model (per-covariate effect curves with posterior
#'   bands), `"tree"` onto a CART partition.
#' @param covariates Optional character or integer vector selecting the
#'   covariates to project onto. Default `NULL` uses all of them. No selection
#'   is performed: choosing this subset is the user's job.
#' @param penalty Regularisation for the linear projection: `"none"` for
#'   ordinary least squares, or `"ridge"`, `"lasso"`, `"elastic_net"`. The
#'   penalised options are the ones that stay defined when the covariate set is
#'   larger than `n`, which is the regime this package targets.
#'   `family = "linear"` only.
#' @param alpha Elastic-net mixing parameter, `penalty = "elastic_net"`.
#'   `alpha = 1` is the lasso, `alpha = 0` the ridge.
#' @param lambda Penalty level for the penalised projections. `NULL` or
#'   `"1se"` cross-validates on the posterior mean and takes `lambda.1se`,
#'   `"min"` takes `lambda.min`, and a number is used directly.
#' @param df_spline Degrees of freedom per natural spline,
#'   `family = "additive"`.
#' @param max_leaves Maximum number of leaves, `family = "tree"`.
#' @param weights Optional non-negative weights over the observations,
#'   defining the population the projection is taken over. Default uniform.
#' @param level Credible level for all reported intervals.
#' @param scale_label Provenance label carried into the printed output
#'   (default method).
#' @param ... Passed on to the default method.
#' @rdname posterior_projection
#'
#' @details
#' The penalised projections are solved **separately for every draw**, at a
#' common `lambda`: for each draw,
#' \eqn{\gamma^{(s)} = \arg\min \|f^{(s)} - X\gamma\|^2 + \lambda\|\gamma\|_1}.
#' The lasso is a nonlinear operator, so the projection of the posterior mean is
#' not the mean of the projections; solving per draw is what gives a genuine
#' posterior over sparse projections, with coefficients that are exactly zero in
#' some draws and not in others. It costs one glmnet solve per draw.
#'
#' `lambda` indexes the summary family, so it is chosen once for all draws, by
#' cross-validating on the posterior mean. Holding out observations asks how
#' well a sparse linear summary generalises across the covariate space, which is
#' a real question even though a posterior draw carries no observation noise.
#'
#' With `penalty = "ridge"` nothing is set to zero, so the result is a
#' regularised projection rather than a lower-dimensional summary.
#'
#' @return An object of class `"PosteriorProjection"`, a list with
#'   `coefficients` (`S` x `q` draws of the projection coefficients),
#'   `coef_summary`, `rho2` (an `S`-vector of summary \eqn{R^2}),
#'   `covariates`, `structure`, and provenance fields.
#' @export
posterior_projection.default <- function(object, X,
                              family       = c("linear", "additive", "tree"),
                              covariates   = NULL,
                              penalty      = c("none", "ridge", "lasso",
                                               "elastic_net"),
                              alpha        = 0.5,
                              lambda       = NULL,
                              df_spline    = 4,
                              max_leaves   = 6,
                              weights      = NULL,
                              level        = 0.95,
                              scale_label  = "response",
                              target       = "f",
                              ...) {

  family  <- match.arg(family)
  penalty <- match.arg(penalty)
  if (penalty != "none" && family != "linear")
    stop("penalty = '", penalty, "' applies only to family = 'linear'.")

  draws <- as.matrix(object)
  if (!all(is.finite(draws)))
    stop("The posterior draws contain ", sum(!is.finite(draws)),
         " non-finite values; the projection would silently propagate them. ",
         "Inspect the stored posterior sample on the fitted object ",
         "(train_predictions_sample / _treat / _control).", call. = FALSE)
  X     <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  n <- nrow(X); S <- nrow(draws)
  if (ncol(draws) != n)
    stop("draws must be S x n with n = nrow(X); got ", ncol(draws), " vs ",
         n, ".")

  if (!is.null(covariates)) {
    idx <- if (is.character(covariates)) match(covariates, colnames(X)) else
      as.integer(covariates)
    if (anyNA(idx))
      stop("These `covariates` are not columns of the design: ",
           paste(covariates[is.na(idx)], collapse = ", "), ".")
    X <- X[, idx, drop = FALSE]
  }
  p_used <- ncol(X)

  w <- if (is.null(weights)) rep(1, n) else as.numeric(weights)
  if (length(w) != n) stop("weights must have length nrow(X).")

  Ft <- t(draws)                                   # n x S
  a  <- (1 - level) / 2

  if (family == "tree") {
    st    <- .pp_tree_structure(rowMeans(Ft), X, w, max_leaves)
    pj    <- .pp_wls(st$B, Ft, w)
    covs  <- st$covariates
    edf   <- ncol(st$B)
    coefs <- t(pj$coef)
    colnames(coefs) <- colnames(st$B)

  } else if (penalty != "none") {
    al <- switch(penalty, ridge = 0, lasso = 1, elastic_net = alpha)
    cv <- .pp_cv_lambda(rowMeans(Ft), X, w, al, lambda)
    pj <- .pp_penalised_draws(Ft, X, w, al, cv$lambda)

    ## Keep Ft aligned with pj$fitted if failed solves were dropped.
    if (length(pj$kept) < ncol(Ft)) {
      Ft <- Ft[, pj$kept, drop = FALSE]
      S  <- ncol(Ft)
    }

    st    <- list(lambda = cv$lambda, rule = cv$rule, cv = cv$cv,
                  nonzero = pj$nonzero)
    covs  <- setdiff(rownames(pj$coef), "(Intercept)")
    edf   <- pj$edf
    coefs <- t(pj$coef)

  } else if (family == "linear") {
    B <- .pp_basis_linear(X)
    .pp_check_dim(ncol(B), n, family, p_used)
    pj    <- .pp_wls(B, Ft, w)
    st    <- list(B = B)
    covs  <- colnames(X)
    edf   <- ncol(B)
    coefs <- t(pj$coef)
    colnames(coefs) <- colnames(B)

  } else {
    ad <- .pp_basis_additive(X, df_spline)
    .pp_check_dim(ncol(ad$B), n, family, p_used)
    pj    <- .pp_wls(ad$B, Ft, w)
    st    <- ad
    covs  <- colnames(X)
    edf   <- ncol(ad$B)
    coefs <- t(pj$coef)
    colnames(coefs) <- colnames(ad$B)
  }

  resid  <- rowMeans(Ft) - rowMeans(pj$fitted)
  df_res <- n - edf

  out <- list(
    family       = family,
    penalty      = penalty,
    alpha        = if (penalty == "elastic_net") alpha else NA_real_,
    covariates   = covs,
    structure    = st,
    coefficients = coefs,
    coef_summary = .pp_coef_summary(coefs, a),
    rho2         = .pp_rho2(Ft, pj$fitted, w),
    residuals    = resid,
    sigma        = sqrt(sum(w * resid^2) / max(df_res, 1)),
    edf          = edf,
    df_residual  = df_res,
    scale        = scale_label,
    target       = target,
    level        = level,
    n_draws      = S,
    n_obs        = n,
    p            = p_used,
    X            = X,
    weights      = w
  )
  class(out) <- "PosteriorProjection"
  out
}


#' Print a posterior projection
#'
#' Laid out like [summary.lm()]: the residual distribution, the coefficient
#' table, and the residual scale and \eqn{R^2}. Uncertainty is reported as
#' posterior credible intervals throughout.
#'
#' @param x A `PosteriorProjection` object.
#' @param digits Significant digits for the coefficient table.
#' @param n_max Maximum number of coefficients to display.
#' @param ... Currently unused.
#' @return Invisibly returns `x`.
#' @export
print.PosteriorProjection <- function(x, digits = 4, n_max = 25, ...) {

  cat("\nProjection residuals:\n")
  res <- x$residuals[is.finite(x$residuals)]
  if (length(res) < length(x$residuals))
    cat("  [", length(x$residuals) - length(res),
        " non-finite residuals omitted]\n", sep = "")
  if (length(res)) {
    rq <- stats::quantile(res)
    names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
    print(signif(rq, digits))
  }

  a  <- (1 - x$level) / 2
  cs <- x$coef_summary
  lab <- paste0(format(100 * c(a, 1 - a), trim = TRUE), "%")

  if (!is.null(x$structure$cv)) {
    cv <- x$structure$cv
    i  <- c(min = which(cv$lambda == cv$lambda.min)[1],
            `1se` = which(cv$lambda == cv$lambda.1se)[1])
    cat("\nMeasure: ", cv$name, "\n", sep = "")
    print(data.frame(Lambda  = signif(cv$lambda[i], 5),
                     Index   = as.integer(i),
                     Measure = signif(cv$cvm[i], 4),
                     SE      = signif(cv$cvsd[i], 4),
                     Nonzero = as.integer(cv$nzero[i]),
                     row.names = c("min", "1se")))
    cat("Using lambda.", x$structure$rule, "\n", sep = "")
  } else if (x$penalty != "none") {
    cat("\nLambda: ", signif(x$structure$lambda, 5), " (supplied)\n", sep = "")
  }

  if (x$family == "additive") {
    ef <- .pp_additive_effects(x, a)
    cat("\nPartial effects (", x$scale,
        "), 10th to 90th percentile of each covariate:\n", sep = "")
    em <- cbind(ef$mean, ef$sd, ef$lower, ef$upper)
    dimnames(em) <- list(ef$variable,
                         c("Estimate", "Post.SD", lab[1], lab[2]))
    print(signif(em, digits))
    cat("Basis coefficients are in x$coef_summary; ",
        "plot(type = \"effects\") draws the curves.\n", sep = "")

  } else if (nrow(cs) <= n_max) {
    cat("\nCoefficients (", x$scale, "):\n", sep = "")
    cm <- cbind(cs$mean, cs$sd, cs$lower, cs$upper)
    nms <- c("Estimate", "Post.SD", lab[1], lab[2])
    if (!is.null(x$structure$nonzero)) {
      cm  <- cbind(cm, x$structure$nonzero)
      nms <- c(nms, "Nonzero")
    }
    dimnames(cm) <- list(cs$variable, nms)
    print(signif(cm, digits))
  } else {
    cat("\nCoefficients: ", nrow(cs), " selected in at least one draw, ",
        "not shown. See x$coef_summary.\n", sep = "")
  }

  r2 <- x$rho2[is.finite(x$rho2)]
  ci <- if (length(r2)) stats::quantile(r2, c(a, 1 - a)) else c(NA, NA)
  cat("\nResidual standard error: ", signif(x$sigma, digits), " on ",
      round(x$df_residual, 1), " degrees of freedom\n", sep = "")
  cat("Summary R-squared: ", signif(mean(r2), digits),
      "  (", round(100 * x$level), "% CI: ", signif(ci[1], digits), ", ",
      signif(ci[2], digits), ")\n", sep = "")
  cat("Draws: ", x$n_draws, "   Observations: ", x$n_obs,
      "   Covariates: ", x$p, "\n", sep = "")

  if (x$family == "tree")
    cat("Leaves: ", x$structure$size, "   Split on: ",
        paste(x$covariates, collapse = ", "), "\n", sep = "")

  invisible(x)
}

#' Lay out an rpart tree for plotting
#'
#' Computes node coordinates for a fitted `rpart` tree so it can be drawn in
#' ggplot2 rather than through `plot.rpart()`'s base graphics. Leaves are
#' placed left to right in tree order, internal nodes at the midpoint of their
#' children, and depth is read off the node id (`floor(log2(id))`).
#'
#' @return A data frame with one row per node: `row` (row of `pr$frame`), `id`
#'   (node id), `x`, `y`, and `leaf`.
#' @noRd
.pp_tree_layout <- function(pr) {
  frame   <- pr$frame
  ids     <- as.integer(rownames(frame))
  row_of  <- stats::setNames(seq_along(ids), as.character(ids))
  isleaf  <- frame$var == "<leaf>"
  xs      <- numeric(length(ids))
  counter <- 0

  assign_x <- function(id) {
    r <- row_of[[as.character(id)]]
    if (isleaf[r]) {
      counter <<- counter + 1
      xs[r]   <<- counter
    } else {
      lx <- assign_x(2L * id)
      rx <- assign_x(2L * id + 1L)
      xs[r] <<- (lx + rx) / 2
    }
    xs[r]
  }
  assign_x(1L)

  data.frame(row = seq_along(ids), id = ids, x = xs,
             y = -floor(log2(ids)), leaf = isleaf)
}

#' Plot a posterior projection
#'
#' @param x A `PosteriorProjection` object.
#' @param type `"coefficients"` is a caterpillar plot of the projection
#'   coefficients; `"rho2"` the posterior of the summary \eqn{R^2}; `"path"`
#'   the summary \eqn{R^2} against \eqn{\log\lambda} (`penalty = "ridge"`
#'   only); `"effects"` the per-covariate partial effect curves with credible
#'   bands (`family = "additive"` only); `"tree"` the CART partition, with the
#'   posterior interval for each leaf (`family = "tree"` only).
#' @param n_max Maximum number of coefficients to show, `type =
#'   "coefficients"`. The largest by absolute posterior mean are kept.
#' @param ... Currently unused.
#' @return A \pkg{ggplot2} object.
#' @export
plot.PosteriorProjection <- function(x, type = c("coefficients", "rho2",
                                                 "path", "effects", "tree"),
                                     n_max = 25, ...) {
  type <- match.arg(type)
  .check_ggplot2()

  if (type == "path") {
    if (is.null(x$structure$cv))
      stop("type = 'path' needs a cross-validated penalty; supply ",
           "lambda = NULL, \"1se\" or \"min\".")
    cv <- x$structure$cv
    df <- data.frame(loglambda = log(cv$lambda), cvm = cv$cvm,
                     lo = cv$cvlo, hi = cv$cvup, nzero = cv$nzero)

    return(
      ggplot2::ggplot(df, ggplot2::aes(x = .data$loglambda, y = .data$cvm)) +
        ggplot2::geom_linerange(
          ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
          colour = "grey70"
        ) +
        ggplot2::geom_point(colour = "firebrick", size = 1.2) +
        ggplot2::geom_vline(xintercept = log(c(cv$lambda.min, cv$lambda.1se)),
                            linetype = "dashed", colour = "grey40") +
        ggplot2::labs(
          x = expression(log(lambda)), y = cv$name,
          title = paste0("Cross-validation (", x$penalty, ")"),
          subtitle = "dashed lines mark lambda.min and lambda.1se"
        ) +
        ggplot2::theme_bw()
    )
  }

  if (type == "tree") {
    if (x$family != "tree") stop("type = 'tree' needs family = 'tree'.")
    pr  <- x$structure$tree
    lay <- .pp_tree_layout(pr)

    kid  <- lay[lay$id != 1L, , drop = FALSE]
    par_ <- lay[match(kid$id %/% 2L, lay$id), , drop = FALSE]
    edges <- data.frame(x = par_$x, y = par_$y, xend = kid$x, yend = kid$y)

    kid$branch <- as.character(labels(pr))[kid$row]
    kid$mx <- (kid$x + par_$x) / 2
    kid$my <- (kid$y + par_$y) / 2

    inner <- lay[!lay$leaf, , drop = FALSE]
    inner$label <- as.character(pr$frame$var[inner$row])

    cs   <- x$coef_summary
    leaf <- lay[match(x$structure$leaf_rows, lay$row), , drop = FALSE]
    leaf$label <- sprintf("%.2f\n[%.2f, %.2f]\nn = %d",
                          cs$mean, cs$lower, cs$upper,
                          pr$frame$n[leaf$row])

    ggplot2::ggplot() +
      ggplot2::geom_segment(
        data = edges,
        ggplot2::aes(x = .data$x, y = .data$y,
                     xend = .data$xend, yend = .data$yend),
        colour = "grey60"
      ) +
      ggplot2::geom_label(
        data = kid,
        ggplot2::aes(x = .data$mx, y = .data$my, label = .data$branch),
        size = 2.4, colour = "grey30", label.size = 0, fill = "white"
      ) +
      ggplot2::geom_label(
        data = inner,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        size = 3, fontface = "bold", fill = "grey95"
      ) +
      ggplot2::geom_label(
        data = leaf,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        size = 2.6, fill = "steelblue", alpha = 0.25, lineheight = 0.95
      ) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.12)) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.08)) +
      ggplot2::labs(
        title = "Projected tree",
        subtitle = paste0("Leaves show the posterior mean and ",
                          round(100 * x$level), "% interval (", x$scale, ")")
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  } else if (type == "rho2") {
    ci <- quantile(x$rho2, c((1 - x$level) / 2, 1 - (1 - x$level) / 2))

    ggplot2::ggplot(data.frame(rho2 = x$rho2),
                    ggplot2::aes(x = .data$rho2)) +
      ggplot2::geom_density(fill = "steelblue", alpha = 0.3,
                            colour = "steelblue") +
      ggplot2::geom_vline(xintercept = median(x$rho2), colour = "steelblue",
                          linewidth = 0.8) +
      ggplot2::geom_vline(xintercept = ci, colour = "steelblue",
                          linetype = "dashed", linewidth = 0.5) +
      ggplot2::labs(
        x = expression(rho^2), y = "Density",
        title = "How much of the model the projection captures",
        subtitle = paste0(round(100 * x$level), "% credible interval")
      ) +
      ggplot2::theme_bw()

  } else if (type == "coefficients") {
    cs <- x$coef_summary
    cs <- cs[cs$variable != "(Intercept)", , drop = FALSE]
    dropped <- max(0L, nrow(cs) - n_max)
    if (dropped > 0)
      cs <- cs[order(abs(cs$mean), decreasing = TRUE)[seq_len(n_max)], ]
    cs$variable <- factor(cs$variable, levels = cs$variable[order(cs$mean)])

    ggplot2::ggplot(cs, ggplot2::aes(x = .data$mean, y = .data$variable)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                          colour = "grey50") +
      ggplot2::geom_linerange(
        ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
        colour = "steelblue", linewidth = 1
      ) +
      ggplot2::geom_point(size = 2, colour = "steelblue") +
      ggplot2::labs(
        x = paste0("Projection coefficient (", x$scale, ")"),
        y = NULL,
        title = "Posterior of the projection coefficients",
        subtitle = paste0(round(100 * x$level), "% credible intervals",
                          if (dropped > 0)
                            paste0("; ", dropped, " smaller coefficients ",
                                   "not shown") else "")
      ) +
      ggplot2::theme_bw()

  } else {
    if (x$family != "additive")
      stop("type = 'effects' is only defined for family = 'additive'.")
    ad <- x$structure
    a  <- (1 - x$level) / 2

    df <- do.call(rbind, lapply(names(ad$index), function(nm) {
      xj   <- x$X[, nm]
      grid <- seq(min(xj), max(xj), length.out = 100)
      ob   <- ad$splines[[nm]]
      Bg   <- if (is.null(ob)) matrix(grid, ncol = 1) else predict(ob, grid)
      cf   <- x$coefficients[, ad$index[[nm]], drop = FALSE]   # S x q
      curves <- Bg %*% t(cf)                                   # grid x S
      curves <- sweep(curves, 2, colMeans(curves))    # centre each draw
      qs <- t(apply(curves, 1, quantile, c(a, 0.5, 1 - a)))
      data.frame(variable = nm, x = grid,
                 lower = qs[, 1], median = qs[, 2], upper = qs[, 3],
                 stringsAsFactors = FALSE)
    }))

    rug_df <- do.call(rbind, lapply(names(ad$index), function(nm)
      data.frame(variable = nm, x = x$X[, nm], stringsAsFactors = FALSE)))

    ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$median)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
        fill = "steelblue", alpha = 0.25
      ) +
      ggplot2::geom_line(colour = "steelblue", linewidth = 0.8) +
      ggplot2::geom_rug(data = rug_df, ggplot2::aes(x = .data$x),
                        inherit.aes = FALSE, alpha = 0.2,
                        length = ggplot2::unit(0.02, "npc")) +
      ggplot2::facet_wrap(~ variable, scales = "free_x") +
      ggplot2::labs(
        x = NULL, y = paste0("Partial effect (", x$scale, ")"),
        title = "Projected additive effects",
        subtitle = paste0(round(100 * x$level),
                          "% credible bands; a projection, not a partial ",
                          "dependence integral")
      ) +
      ggplot2::theme_bw()
  }
}

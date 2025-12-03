#' Wrapper for C++ cwbart implementation
#'
#' @export
cwbart <- function(
    x.train, y.train, 
    x.test = matrix(0.0, 0, 0),
    sparse = FALSE, theta = 0, omega = 1,
    a = 0.5, b = 1, augment = FALSE, rho = NULL,
    usequants = FALSE,
    cont = FALSE, rm.const = TRUE,
    sigest = NA, sigdf = 3, sigquant = .90,
    k = 2.0, power = 2.0, base = .95,
    sigmaf = NA, lambda = NA,
    fmean = mean(y.train),
    w = rep(1, length(y.train)),
    ntree = 200L, numcut = 100L,
    ndpost = 1000L, nskip = 100L, keepevery = 1L,
    nkeeptrain = ndpost, nkeeptest = ndpost,
    nkeeptestmean = ndpost, nkeeptreedraws = ndpost,
    printevery = 100L, transposed = FALSE
) {

    ## ------------------------
    ## 1. Prepare data
    ## ------------------------
    n = length(y.train)

    if (!transposed) {
        temp = bartModelMatrix(
            x.train, numcut,
            usequants = usequants,
            cont = cont,  rm.const = rm.const
        )
        x.train = t(temp$X)
        numcut = temp$numcut

        if (length(x.test) > 0) {
            x.test = bartModelMatrix(x.test)
            x.test = t(x.test[, temp$rm.const])
        }
        rm.const <- temp$rm.const
        grp <- temp$grp
        rm(temp)
    } else {
        rm.const <- NULL
        grp <- NULL
    }

    if (n != ncol(x.train))
        stop("Length of y.train must equal number of rows in x.train")

    p  = nrow(x.train)
    np = ncol(x.test)

    if (is.null(rho))
        rho <- p

    y.train <- y.train - fmean

    ## ------------------------
    ## 2. Set sigma prior
    ## ------------------------
    nu = sigdf

    if (is.na(lambda)) {
        if (is.na(sigest)) {
            if (p < n) {
                df <- data.frame(t(x.train), y.train)
                sigest <- summary(lm(y.train ~ ., df))$sigma
            } else {
                sigest <- sd(y.train)
            }
        }
        qchi <- qchisq(1 - sigquant, nu)
        lambda <- (sigest^2 * qchi) / nu
    } else {
        sigest <- sqrt(lambda)
    }

    ## tau (prior on mu)
    if (is.na(sigmaf)) {
        tau <- (max(y.train) - min(y.train)) / (2 * k * sqrt(ntree))
    } else {
        tau <- sigmaf / sqrt(ntree)
    }

    ## ------------------------
    ## 3. Call the C++ sampler
    ## ------------------------
    res <- .Call(
        `_ShrinkageTrees_cwbart`,
        ntrain = n,
        ip = p,
        inp = np,
        ixSEXP = x.train,
        iySEXP = y.train,
        ixpSEXP = x.test,
        im = ntree,
        inc = numcut,
        ind = ndpost * keepevery,
        iburn = nskip,
        ipower = power,
        ibase = base,
        itau = tau,
        inu = nu,
        ilambda = lambda,
        isigest = sigest,
        iwSEXP = w,
        idart = sparse,
        itheta = theta,
        ia = a,
        ib = b,
        irho = rho,
        iaug = augment,
        inkeeptrain = nkeeptrain,
        inkeeptest = nkeeptest,
        inkeeptestme = nkeeptestmean,
        inkeeptreedraws = nkeeptreedraws,
        inprintevery = printevery
    )

    ## ------------------------
    ## 4. Post-processing
    ## ------------------------
    res$mu <- fmean
    res$yhat.train.mean <- res$yhat.train.mean + fmean
    res$yhat.train       <- res$yhat.train + fmean
    res$yhat.test.mean   <- res$yhat.test.mean + fmean
    res$yhat.test        <- res$yhat.test + fmean

    attr(res, "class") <- "wbart"

    return(res)
}

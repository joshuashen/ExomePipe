# R functions for generalized linear modeling with independent normal, t, or
# Cauchy prior distribution for the coefficients

# Default prior distribution is Cauchy with center 0 and scale 2.5 for all
# coefficients (except for the intercept, which has a prior scale of 10),
# as described in the paper,
# "A default prior distribution for logistic and other regression models,"
# by Andrew Gelman, Aleks Jakulin, Maria Grazia Pittau, and Yu-Sung Su

# R functions by Andrew Gelman and Yu-Sung Su

# Update 12 Sep 2006

# Run it just like glm() with these additional arguments:
#   prior.mean (default is 0)
#   prior.scale (default is 10 for the intercept and 2.5 for all other coefs)
#   prior.df for t distribution (default is 1 (Cauchy))
# The program is a simple alteration of glm() that uses an approximate EM
# algorithm to update the betas at each step using an augmented regression
# to represent the prior information

# It is recommended that continuous predictors be rescaled by dividing by
# two standard deviations (see the function standardize.R at
# http://www.stat.columbia.edu/~gelman/standardize/standardize.R)

# Examples appear after the functions bayesglm() and bayesglm.fit()

### downloaded from: http://www.stat.columbia.edu/~gelman/standardize/bayesglm.R 

bayesglm <-
function (formula, family = gaussian, data, weights, subset, 
    na.action, start = NULL, etastart, mustart, offset, control = glm.control(...), 
    model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL,
    prior.mean=0, prior.scale=2.5, prior.df=1, n.iter=50, ...) 
{
    call <- match.call()
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "etastart", "mustart", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    switch(method, model.frame = return(mf), glm.fit = 1, stop("invalid 'method': ", 
        method))
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) 
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0)
    weights <- model.weights(mf)
    offset <- model.offset(mf)
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    if (!is.null(offset) && length(offset) != NROW(Y)) 
        stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
            length(offset), NROW(Y)), domain = NA)
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
# If only one set of prior parameters is set, expand to all the coefs
    if (length(prior.mean)==1) prior.mean <- rep(prior.mean, NCOL(X))
    if (length(prior.scale)==1) {
      prior.scale <- rep(prior.scale, NCOL(X))
      if (attr(mt, "intercept") > 0) prior.scale[1] <- 10
    }
    if (length(prior.df)==1) prior.df <- rep(prior.df, NCOL(X))
# change glm.fit to bayesglm.fit
    fit <- bayesglm.fit(x = X, y = Y, weights = weights, start = start, 
        etastart = etastart, mustart = mustart, offset = offset,                    
        family = family, control=glm.control(maxit=n.iter),
        intercept = attr(mt, "intercept") > 0,          
        prior.mean=prior.mean, prior.scale=prior.scale, prior.df=prior.df)          
    if (any(offset) && attr(mt, "intercept") > 0) {
# change glm.fit to bayesglm.fit
      cat ("bayesglm not yet set up to do deviance comparion here\n")
        fit$null.deviance <- bayesglm.fit(x = X[, "(Intercept)", drop = FALSE], 
            y = Y, weights = weights, offset = offset, family = family, 
            control = control, intercept = TRUE,
            prior.mean=prior.mean, prior.scale=prior.scale, prior.df=prior.df)$deviance
    }
    if (model) 
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x) 
        fit$x <- X
    if (!y) 
        fit$y <- NULL
    fit <- c(fit, list(call = call, formula = formula, terms = mt, 
        data = data, offset = offset, control = control, method = method, 
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
            mf)))
    class(fit) <- c("glm", "lm")
    fit
}

bayesglm.fit <- 
function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL, 
    mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
    control = glm.control(...), intercept = TRUE,
          prior.mean=rep(0,NCOL(x)), prior.scale=rep(Inf,NCOL(x)),
          prior.df=rep(Inf,NCOL(x))) 
{
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2]]
    ynames <- if (is.matrix(y)) 
        rownames(y)
    else names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights)) 
        weights <- rep.int(1, nobs)
    if (is.null(offset)) 
        offset <- rep.int(0, nobs)
    variance <- family$variance
    dev.resids <- family$dev.resids
    aic <- family$aic
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    if (!is.function(variance) || !is.function(linkinv)) 
        stop("'family' argument seems not to be a valid family object")
    valideta <- family$valideta
    if (is.null(valideta)) 
        valideta <- function(eta) TRUE
    validmu <- family$validmu
    if (is.null(validmu)) 
        validmu <- function(mu) TRUE
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta)) 
            stop("invalid linear predictor values in empty model")
        mu <- linkinv(eta)
        if (!validmu(mu)) 
            stop("invalid fitted means in empty model")
        dev <- sum(dev.resids(y, mu, weights))
        w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric(0)
        iter <- 0
    }
    else {
        coefold <- NULL
        eta <- if (!is.null(etastart)) 
            etastart
        else if (!is.null(start)) 
            if (length(start) != nvars) 
                stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s", 
                  nvars, paste(deparse(xnames), collapse = ", ")), 
                  domain = NA)
            else {
                coefold <- start
                offset + as.vector(if (NCOL(x) == 1) 
                  x * start
                else x %*% start)
            }
        else family$linkfun(mustart)
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta))) 
            stop("cannot find valid starting values: please specify some")
        devold <- sum(dev.resids(y, mu, weights))
        boundary <- conv <- FALSE
# Initialize the prior sd at the prior scale
        prior.sd <- prior.scale
#
        for (iter in 1:control$maxit) {
            good <- weights > 0
            varmu <- variance(mu)[good]
            if (any(is.na(varmu))) 
                stop("NAs in V(mu)")
            if (any(varmu == 0)) 
                stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            if (any(is.na(mu.eta.val[good]))) 
                stop("NAs in d(mu)/d(eta)")
            good <- (weights > 0) & (mu.eta.val != 0)
            if (all(!good)) {
                conv <- FALSE
                warning("no observations informative at iteration ", 
                  iter)
                break
            }
            z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
            w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
            ngoodobs <- as.integer(nobs - sum(!good))
# This is where we augment the data with the prior information
            z.star <- c (z, prior.mean)
            x.star <- rbind (x, diag(NCOL(x)))
            w.star <- c (w, 1/prior.sd)
            good.star <- c (good, rep(TRUE,NCOL(x)))
            ngoodobs.star <- ngoodobs + NCOL(x)
            fit <- .Fortran("dqrls", qr = x.star[good.star, ] * w.star, n = ngoodobs.star, 
                p = nvars, y = w.star * z.star, ny = as.integer(1), tol = min(1e-07, 
                  control$epsilon/1000), coefficients = double(nvars), 
                residuals = double(ngoodobs.star), effects = double(ngoodobs.star), 
                rank = integer(1), pivot = 1:nvars, qraux = double(nvars), 
                work = double(2 * nvars), PACKAGE = "base")
            if (any(!is.finite(fit$coefficients))) {
                conv <- FALSE
                warning("non-finite coefficients at iteration ", 
                  iter)
                break
            }
# Now update the prior scale
            prior.sd <- ifelse (prior.df==Inf, prior.scale, sqrt (
              (fit$coefficients^2 + prior.df*prior.scale^2)/(1 + prior.df)))
# Comment out the next 3 lines because all is ok with a Bayesian prior dist.
#            if (nobs < fit$rank) 
#                stop(gettextf("X matrix has rank %d, but only %d observations", 
#                  fit$rank, nobs), domain = NA)
            start[fit$pivot] <- fit$coefficients
            eta <- drop(x %*% start)
            mu <- linkinv(eta <- eta + offset)
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace) 
                cat("Deviance =", dev, "Iterations -", iter, 
                  "\n")
            boundary <- FALSE
            if (!is.finite(dev)) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found: please supply starting values", 
                    call. = FALSE)
                warning("step size truncated due to divergence", 
                  call. = FALSE)
                ii <- 1
                while (!is.finite(dev)) {
                  if (ii > control$maxit) 
                    stop("inner loop 1; cannot correct step size")
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                  dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }
            if (!(valideta(eta) && validmu(mu))) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found: please supply starting values", 
                    call. = FALSE)
                warning("step size truncated: out of bounds", 
                  call. = FALSE)
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                  if (ii > control$maxit) 
                    stop("inner loop 2; cannot correct step size")
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                }
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }
            if (iter>1 & abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
                conv <- TRUE
                coef <- start
                break
            }
            else {
                devold <- dev
                coef <- coefold <- start
            }
        }
        if (!conv) 
            warning("algorithm did not converge")
        if (boundary) 
            warning("algorithm stopped at boundary value")
        eps <- 10 * .Machine$double.eps
        if (family$family == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps)) 
                warning("fitted probabilities numerically 0 or 1 occurred")
        }
        if (family$family == "poisson") {
            if (any(mu < eps)) 
                warning("fitted rates numerically 0 occurred")
        }
        if (fit$rank < nvars) 
            coef[fit$pivot][seq(fit$rank + 1, nvars)] <- NA
        xxnames <- xnames[fit$pivot]
        residuals <- rep.int(NA, nobs)
        residuals[good] <- z - (eta - offset)[good]
        fit$qr <- as.matrix(fit$qr)
        nr <- min(sum(good), nvars)
        if (nr < nvars) {
            Rmat <- diag(nvars)
            Rmat[1:nr, 1:nvars] <- fit$qr[1:nr, 1:nvars]
        }
        else Rmat <- fit$qr[1:nvars, 1:nvars]
        Rmat <- as.matrix(Rmat)
        Rmat[row(Rmat) > col(Rmat)] <- 0
        names(coef) <- xnames
        colnames(fit$qr) <- xxnames
        dimnames(Rmat) <- list(xxnames, xxnames)
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    #if (!EMPTY) 
    #    names(fit$effects) <- c(xxnames[seq(len = fit$rank)], 
    #        rep.int("", sum(good) - fit$rank))
    wtdmu <- if (intercept) 
        sum(weights * y)/sum(weights)
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- if (EMPTY) 
        0
    else fit$rank
    resdf <- n.ok - rank
    aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
    list(coefficients = coef, residuals = residuals, fitted.values = mu, 
        effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat, 
        rank = rank, qr = if (!EMPTY) structure(fit[c("qr", "rank", 
            "qraux", "pivot", "tol")], class = "qr"), family = family, 
        linear.predictors = eta, deviance = dev, aic = aic.model, 
        null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights, 
        df.residual = resdf, df.null = nulldf, y = y, converged = conv, 
        boundary = boundary)
}


# Now, for some examples.  First make sure you have the functions at http://www.stat.columbia.edu/~gelman/bugsR/regression.R

if (0){
  n <- 100
  x1 <- rnorm (n)
  x1 <- (x1-mean(x1))/(2*sd(x1))   # standardization
  x2 <- rbinom (n, 1, .5)
  b0 <- 1
  b1 <- 1.5
  b2 <- 2
  y <- rbinom (n, 1, invlogit(b0+b1*x1+b2*x2))

  M1 <- glm (y ~ x1 + x2, family=binomial(link="logit"))
  display (M1)  # (using the display() function from regression.R)

  M2 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"), prior.scale=Inf, prior.df=Inf)
  display (M2)  # just a test:  this should be identical to classical logit

  M3 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"))  # default Cauchy prior with scale 2.5
  display (M3)

  M4 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"), prior.scale=2.5, prior.df=1)  # Same as M3, explicitly specifying Cauchy prior with scale 2.5
  display (M4)

  M5 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"), prior.scale=2.5, prior.df=7)   # t_7 prior with scale 2.5
  display (M5)

  M6 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"), prior.scale=2.5, prior.df=Inf)  # normal prior with scale 2.5
  display (M6)

# Create separation:  set y=1 whenever x2=1
# Now it should blow up without the prior!

  y <- ifelse (x2==1, 1, y)

  M1 <- glm (y ~ x1 + x2, family=binomial(link="logit"))
  display (M1)

  M2 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"), prior.scale=Inf, prior.df=Inf) # Same as M1
  display (M2)

  M3 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"))
  display (M3)

  M4 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"), prior.scale=2.5, prior.df=1)  # Same as M3
  display (M4)

  M5 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"), prior.scale=2.5, prior.df=7)
  display (M5)

  M6 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"), prior.scale=2.5, prior.df=Inf)
  display (M6)
}

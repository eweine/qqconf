#' QQ Plot with Simultaneous and Pointwise Testing Bounds.
#'
#' Create a qq-plot with with a shaded simultaneous acceptance region and,
#' optionally, lines for a point-wise region. The observed values are plotted
#' against their expected values had they come from the specified distribution.
#' 
#' If any of the points of the qq-plot fall outside the simultaneous acceptance region for the selected
#' level alpha test, that means that we can reject the null hypothesis that the data are i.i.d. draws from the
#' specified distribution. If 'difference' is set to TRUE, the vertical axis plots the 
#' observed quantile minus expected quantile. Set pw.lty to a non-zero line type to plot
#' the pointwise bounds. If pointwise bands are used, then on average, alpha * n of the points will fall outside
#' the bounds under the null hypothesis, so the chance that the qq-plot has any points falling outside of the pointwise bounds
#' is typically much higher than alpha under the null hypothesis. For this reason, a simultaneous region is preferred. 
#' 
#' @param obs The observed data.
#' @param distribution The quantile function for the specified distribution. Defaults to qnorm.
#' Custom distributions are allowed so long as all parameters are supplied in dparams.
#' @param method Method for simultaneous testing bands. Must be either "ell", which applies a level \eqn{\eta} pointwise
#' test to each order statistic such that the Type I error of the global test is \eqn{\alpha}, or "ks" to apply a 
#' Kolmogorov-Smirnov test. For \eqn{\alpha} = .01, .05, and .1, "ell" is recommended.
#' @param alpha Type I error of global test of if the data comes from the reference distribution.
#' @param difference Whether to plot the difference between the observed and
#'   expected values on the vertical axis.
#' @param log10 Whether to plot axes on -log10 scale (e.g. to see small p-values). Can only be used for strictly
#' positive distributions.
#' @param shade.col What color to use for the simultaneous acceptance region.
#' @param add Whether to add points to an existing plot. 
#' @param dparams List of additional parameters for the quantile function of the distribution
#'   (e.g. df=1). Will be estimated if not provided and an appropriate estimation procedure exists.
#'   For the normal distribution, we estimate the mean as the median and the standard deviation as \eqn{Sn} from the paper by Rousseeuw and Croux 1993
#'   "Alternatives to the Median Absolute Deviation". For all other distributions,
#'   the code uses MLE to estimate the parameters. Note that estimation is not implemented for custom distributions, so all
#'   parameters of the distribution must be provided by the user.
#' @param bounds_params List of optional parameters for get_bounds_two_sided
#'   (i.e. tol, max_it, method).
#' @param pw.lty Line type for the pointwise error bounds. Set to non-zero for a line.
#' @param pw.col Color for the pointwise bounds line.
#' @param ... Additional parameters for the plot.
#' 
#' @export
#'
#' @examples
#' x <- rchisq(1000, 1)
#' qq_conf_plot(x, qchisq, dparams=list(df=1), pw.lty=3) # Plots x against a 1-df chisquare
#'
#' y <- runif(893)
#' qq_conf_plot(y, difference = TRUE, log10 = TRUE, bounds_params = list(method = "search"), pch=3)
qq_conf_plot <- function(obs,
                         distribution = qnorm,
                         method = c("ell", "ks"),
                         alpha = 0.05,
                         difference = FALSE,
                         log10 = FALSE,
                         shade.col = 'gray',
                         add = FALSE,
                         dparams = NULL,
                         bounds_params = NULL,
                         pw.lty = 0,
                         pw.col = 'black',
                         ...) {
  
  dist_name <- as.character(substitute(distribution))
  
  if(is.null(dparams)) {
    # equivalence between base R and MASS::fitdistr distribution names
    corresp <- function(distributionName) {
      switch(
        distributionName,
        qbeta = "beta",
        qcauchy = "cauchy",
        qchisq = "chi-squared",
        qexp = "exponential",
        qf = "f",
        qgamma = "gamma",
        qgeom = "geometric",
        qlnorm = "log-normal",
        qlogis = "logistic",
        qnorm = "normal",
        qnbinom = "negative binomial",
        qpois = "poisson",
        qt = "t",
        qweibull = "weibull",
        NULL
      )
    }
    
    # initial value for some distributions
    initVal <- function(distributionName) {
      switch(
        distributionName,
        qbeta = list(shape1 = 1, shape2 = 1),
        qchisq = list(df = 1),
        qf = list(df1 = 1, df2 = 2),
        qt = list(df = 1),
        NULL
      )
    }
    
    suppressWarnings({
      if(!is.null(corresp(dist_name))) {
        if(is.null(initVal(dist_name))) {
          if(corresp(dist_name) == "normal") {
            
            # Use special estimators for the normal distribution
            dparams <- c()
            dparams['mean'] <- median(x = obs)
            dparams['sd'] <- robustbase::Sn(x = obs)
            
          }
          dparams <- MASS::fitdistr(x = obs, densfun = corresp(dist_name))$estimate
        } else {
          dparams <- MASS::fitdistr(x = obs, densfun = corresp(dist_name), start = initVal(dist_name))$estimate
        }
      }
    })
  }
  
  dots <- list(...)
  method <- match.arg(method)
  if ( is.null(dots$ylab)) {
    if (difference && log10) {
      ylab <- expression("-log"[10]*"(Observed quantiles) + log"[10]*"(Expected quantiles)")
    } else if(difference) {
      ylab <- 'Obseved quantiles - Expected quantiles'
    } else if (log10) {
      ylab <- expression("-log"[10]*"(Observed quantiles)")
    } else {
      ylab <- 'Observed quantiles'
    }
  } else {
    ylab <- dots$ylab
    dots <- within(dots, rm(ylab))
  }
  if (is.null(dots$xlab)) {
    if (log10) {
      xlab <- expression("-log"[10]*"(Expected quantiles)")
    } else {
      xlab <- 'Expected quantiles'
    }
  } else {
    xlab <- dots$xlab
    dots <- within(dots, rm(xlab))
  }
  
  samp.size <- length(obs)
  conf.int <- 1 - alpha
  conf <- c(alpha / 2, conf.int + alpha / 2)
  ## The observed and expected quantiles. Expected quantiles are based on the specified
  ## distribution
  obs.pts <- sort(obs)
  exp.pts <- do.call(distribution, c(list(p=ppoints(samp.size, a=0)), dparams))
  if (log10 == TRUE) {
    exp.pts <- -log10(exp.pts)
    if (any(obs.pts <= 0)) {
      
      stop("log10 scaling can only be used with strictly positive distributions.")
      
    }
    obs.pts <- -log10(obs.pts)
  }
  if (difference) {
    y.pts <- obs.pts - exp.pts
  } else {
    y.pts <- obs.pts
  }
  
  ## When not adding points to a qq-plot compute pointwise and global confidence bounds.
  if (!add) {
    left <- exp.pts[1]
    right <- exp.pts[samp.size]
    bottom <- min(y.pts) #obs.pts[1]
    top <- max(y.pts) #obs.pts[samp.size]
    do.call(plot, c(list(x=c(left, right), y=c(bottom, top), type='n', xlab=xlab, ylab=ylab), dots))
    pointwise.low <- do.call(distribution,
                             c(list(p=qbeta(conf[1], 1:samp.size, samp.size:1)), dparams))
    pointwise.high <- do.call(distribution,
                              c(list(p=qbeta(conf[2], 1:samp.size, samp.size:1)), dparams))
    
    global.bounds <- do.call(get_bounds_two_sided,
                             c(list(alpha = alpha, n = samp.size), bounds_params))
    # Here, have to figure out how to do this for the KS test
    # I don't think that this should be too hard, but I'm not completely sure
    if (method == "ell") {
      
      global.low <- do.call(distribution, c(list(p = global.bounds$lower_bound), dparams))
      global.high <- do.call(distribution, c(list(p = global.bounds$upper_bound), dparams))
      
    } else if (method == "ks") {
      
      probs <- ppoints(samp.size)
      epsilon <- sqrt((1 / (2 * samp.size)) * log(2 / (1 - conf.int)))
      lp <- pmax(probs - epsilon, rep(0, samp.size))
      up <- pmin(probs + epsilon, rep(1, samp.size))
      global.low <- do.call(distribution, c(list(p = lp), dparams))
      global.high <- do.call(distribution, c(list(p = up), dparams))
      
    }

    # code to extend region for visibility
    global.low <- c(global.low[1], global.low, global.low[samp.size])
    global.high <- c(global.high[1], global.high, global.high[samp.size])
    pointwise.low <- c(pointwise.low[1], pointwise.low, pointwise.low[samp.size])
    pointwise.high <- c(pointwise.high[1], pointwise.high, pointwise.high[samp.size])
    c <- .5
    low_exp_pt <- c * do.call(distribution, c(list(p=c(1 / (samp.size + 2))))) + (1 - c) * exp.pts[1]
    high_exp_pt <- c * do.call(distribution, c(list(p=c((samp.size + 1) / (samp.size + 2))))) + (1 - c) * exp.pts[samp.size]
    exp.pts <- c(low_exp_pt, exp.pts, high_exp_pt)
    
    if (log10 == TRUE) {
      pointwise.low <- -log10(pointwise.low)
      pointwise.high <- -log10(pointwise.high)
      global.low <- -log10(global.low)
      global.high <- -log10(global.high)
    }
    if (difference) {
      polygon(c(exp.pts, rev(exp.pts)),
              c(global.low - exp.pts, rev(global.high) - rev(exp.pts)),
              border=NA, col=shade.col)
      lines(exp.pts, pointwise.low - exp.pts, lty = pw.lty, col = pw.col, ...)
      lines(exp.pts, pointwise.high - exp.pts, lty = pw.lty, col = pw.col, ...)
    } else {
      bottom <- min(y.pts) 
      top <- max(y.pts)
      polygon(c(exp.pts, rev(exp.pts)),
              pmin(pmax(c(global.low, rev(global.high)), bottom), top),
              border=NA, col=shade.col)
      lines(exp.pts, pointwise.low, lty = pw.lty, col = pw.col, ...)
      lines(exp.pts, pointwise.high, lty = pw.lty, col = pw.col, ...)
    }
  }
  points(exp.pts[1:samp.size], y.pts, ...)
  if (difference) {
    abline(h = 0, ...)
  } else {
    abline(0, 1, ...)
  }
  
}
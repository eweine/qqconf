#' PP Plot with Simultaneous and Pointwise Testing Bounds.
#'
#' Create a pp-plot with with a shaded simultaneous acceptance region and,
#' optionally, lines for a point-wise region. The observed values are plotted
#' against their expected values had they come from the specified distribution.
#' 
#' If any of the points of the pp-plot fall outside the simultaneous acceptance region for the selected
#' level alpha test, that means that we can reject the null hypothesis that the data are i.i.d. draws from the
#' specified distribution. If 'difference' is set to TRUE, the vertical axis plots the 
#' observed distribution minus expected distribution. Set pw.lty to a non-zero line type to plot
#' the pointwise bounds. If pointwise bands are used, then on average, alpha * n of the points will fall outside
#' the bounds under the null hypothesis, so the chance that the pp-plot has any points falling outside of the pointwise bounds
#' is typically much higher than alpha under the null hypothesis. For this reason, a simultaneous region is preferred. 
#' 
#' @param obs The observed data.
#' @param distribution The distribution function for the specified distribution. Defaults to pnorm.
#' Custom distributions are allowed so long as all parameters are supplied in dparams.
#' @param method Method for simultaneous testing bands. Must be either "ell", which applies a level \eqn{\eta} pointwise
#' test to each order statistic such that the Type I error of the global test is \eqn{\alpha}, or "ks" to apply a 
#' Kolmogorov-Smirnov test. For \eqn{\alpha} = .01, .05, and .1, "ell" is recommended.
#' @param alpha Type I error of global test of if the data comes from the reference distribution.
#' @param difference Whether to plot the difference between the observed and
#'   expected values on the vertical axis.
#' @param log10 Whether to plot axes on -log10 scale (e.g. to see small p-values).
#' @param shade.col What color to use for the simultaneous acceptance region.
#' @param add Whether to add points to an existing plot. 
#' @param dparams List of additional parameters for the distribution function of the distribution
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
#' pp_conf_plot(x, qchisq, dparams=list(df=1), pw.lty=3) # Plots x against a 1-df chisquare
#'
#' y <- runif(893)
#' pp_conf_plot(y, difference = TRUE, log10 = TRUE, bounds_params = list(method = "search"), pch=3)
pp_conf_plot <- function(obs,
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
  
  if(is.null(dparams)) {
    # equivalence between base R and MASS::fitdistributionr distributionribution names
    corresp <- function(distributionName) {
      switch(
        distributionName,
        pbeta = "beta",
        pcauchy = "cauchy",
        pchisq = "chi-squared",
        pexp = "exponential",
        pf = "f",
        pgamma = "gamma",
        pgeom = "geometric",
        plnorm = "log-normal",
        plogis = "logistic",
        pnorm = "normal",
        pnbinom = "negative binomial",
        ppois = "poisson",
        pt = "t",
        pweibull = "weibull",
        NULL
      )
    }
    
    # initial value for some distributionributions
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
      if(!is.null(corresp(distributionribution))) {
        if(is.null(initVal(distributionribution))) {
          if(corresp(distributionribution) == "normal") {
            
            # Use special estimators for the normal distributionribution
            dparams <- c()
            dparams['mean'] <- median(x = smp)
            dparams['sd'] <- robustbase::Sn(x = smp)
            
          }
          dparams <- MASS::fitdistributionr(x = smp, densfun = corresp(distributionribution))$estimate
        } else {
          dparams <- MASS::fitdistributionr(x = smp, densfun = corresp(distributionribution), start = initVal(distributionribution))$estimate
        }
      }
    })
  }
  
  dots <- list(...)
  method <- match.arg(method)
  if ( is.null(dots$ylab)) {
    if (difference && log10) {
      ylab <- expression("-log"[10]*"(Observed distributions) + log"[10]*"(Expected distributions)")
    } else if(difference) {
      ylab <- 'Obseved distributions - Expected distributions'
    } else if (log10) {
      ylab <- expression("-log"[10]*"(Observed distributions)")
    } else {
      ylab <- 'Observed distributions'
    }
  } else {
    ylab <- dots$ylab
    dots <- within(dots, rm(ylab))
  }
  if (is.null(dots$xlab)) {
    if (log10) {
      xlab <- expression("-log"[10]*"(Expected distributions)")
    } else {
      xlab <- 'Expected distributions'
    }
  } else {
    xlab <- dots$xlab
    dots <- within(dots, rm(xlab))
  }
  
  samp.size <- length(obs)
  conf.int <- 1 - alpha
  conf <- c(alpha / 2, conf.int + alpha / 2)
  ## The observed and expected distributions. Expected distributions are based on the specified
  ## distribution
  obs.pts <- sort(obs)
  exp.pts <- do.call(distribution, c(list(p=ppoints(samp.size, a=0)), dparams))
  if (log10 == TRUE) {
    exp.pts <- -log10(exp.pts)
    obs.pts <- -log10(obs.pts)
  }
  if (difference) {
    y.pts <- obs.pts - exp.pts
  } else {
    y.pts <- obs.pts
  }
  
  ## When not adding points to a pp-plot compute pointwise and global confidence bounds.
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
      
      probs <- ppoints(n)
      epsilon <- sqrt((1 / (2 * n)) * log(2 / (1 - conf.int)))
      lp <- pmax(probs - epsilon, rep(0, n))
      up <- pmin(probs + epsilon, rep(1, n))
      lower <- intercept + slope * do.call(qFunc, c(list(p = lp), dparams))
      upper <- intercept + slope * do.call(qFunc, c(list(p = up), dparams))
      
    }
    
    
    if (log10 == TRUE) {
      pointwise.low <- -log10(pointwise.low)
      pointwise.high <- -log10(pointwise.high)
      global.low <- -log10(global.low)
      global.high <- -log10(global.high)
    }
    if (difference) {
      polygon(c(exp.pts, exp.pts[samp.size:1]),
              c(global.low - exp.pts, global.high[samp.size:1] - exp.pts[samp.size:1]),
              border=NA, col=shade.col)
      lines(exp.pts, pointwise.low - exp.pts, lty = pw.lty, col = pw.col, ...)
      lines(exp.pts, pointwise.high - exp.pts, lty = pw.lty, col = pw.col, ...)
    } else {
      polygon(c(exp.pts, exp.pts[samp.size:1]),
              c(global.low, global.high[samp.size:1]),
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
#' PP Plot with Simultaneous and Pointwise Testing Bounds.
#'
#' Create a pp-plot with with a shaded simultaneous acceptance region and,
#' optionally, lines for a point-wise region. The observed values are plotted
#' against their expected values had they come from the specified distribution.
#' 
#' If any of the points of the pp-plot fall outside the simultaneous acceptance region for the selected
#' level alpha test, that means that we can reject the null hypothesis that the data are i.i.d. draws from the
#' specified distribution. If \code{difference} is set to TRUE, the vertical axis plots the 
#' observed probability minus expected probability. If pointwise bounds are used, then on average, alpha * n of the points will fall outside
#' the bounds under the null hypothesis, so the chance that the pp-plot has any points falling outside of the pointwise bounds
#' is typically much higher than alpha under the null hypothesis. For this reason, a simultaneous region is preferred. 
#' 
#' @param obs The observed data.
#' @param distribution The probability function for the specified distribution. Defaults to \code{pnorm}.
#' Custom distributions are allowed as long as all parameters are supplied in dparams.
#' @param method Method for simultaneous testing bands. Must be either "ell" (equal local levels test), which applies a level \eqn{\eta} pointwise
#' test to each order statistic such that the Type I error of the global test is \code{alpha}, or "ks" to apply a 
#' Kolmogorov-Smirnov test. For \code{alpha} = .01, .05, and .1, "ell" is recommended.
#' @param alpha Type I error of global test of whether the data come from the reference distribution.
#' @param difference Whether to plot the difference between the observed and
#'   expected values on the vertical axis.
#' @param log10 Whether to plot axes on -log10 scale (e.g. to see small p-values). 
#' @param right_tail This parameter is only used if \code{log10} is \code{TRUE}. When \code{TRUE},
#' the x-axis is -log10(1 - Expected Probability) and the y-axis is -log10(1 - Observed Probability).
#' When \code{FALSE} (default) the x-axis is -log10(Expected Probability) and the y-axis is 
#' -log10(Observed Probability). The parameter should be set to \code{TRUE} to make
#' observations in the right tail of the distribution easier to see, and set to false to make the 
#' observations in the left tail of the distribution easier to see.
#' @param add Whether to add points to an existing plot. 
#' @param dparams List of additional parameters for the probability function of the distribution
#'   (e.g. df=1). Note that if any parameters of the distribution are specified, parameter estimation will not be performed
#'   on the unspecified parameters, and instead they will take on the default values set by the distribution function. 
#'   For the uniform distribution, parameter estimation is not performed, and
#'   the default parameters are max = 1 and min = 0.
#'   For other distributions parameters will be estimated if not provided.
#'   For the normal distribution, we estimate the mean as the median and the standard deviation as \eqn{Sn} from the paper by Rousseeuw and Croux 1993
#'   "Alternatives to the Median Absolute Deviation". For all other distributions besides uniform and normal,
#'   the code uses MLE to estimate the parameters. Note that estimation is not implemented for custom distributions, so all
#'   parameters of the distribution must be provided by the user.
#' @param bounds_params List of optional parameters for get_bounds_two_sided
#'   (i.e. \code{tol}, \code{max_it}, \code{method}).
#' @param line_params Parameters passed to the line function to modify the line that indicates a perfect fit of the
#'   reference distribution.
#' @param plot_pointwise Boolean indicating whether pointwise bounds should be added to the plot
#' @param pointwise_lines_params Parameters passed to the \code{lines} function that modifies pointwise bounds when plot_pointwise is
#'   set to TRUE.
#' @param points_params Parameters to be passed to the \code{points} function to plot the data.
#' @param polygon_params Parmeters to be passed to the polygon function to construct simultaneous confidence region.
#'   By default \code{border} is set to NA and \code{col} is set to grey.
#' @param ... Additional parameters passed to the plot function.
#' 
#' @export
#' 
#' @return None, PP plot is produced.
#'
#' @examples
#' set.seed(0)
#' smp <- rnorm(100)
#' 
#' # Plot PP plot against normal distribution with mean and variance estimated
#' pp_conf_plot(
#'   obs=smp, 
#'   distribution = pnorm
#' )
#' 
#' # Make same plot on -log10 scale to highlight the left tail,
#' # with radius of plot circles also reduced by .5
#' pp_conf_plot(
#'   obs=smp, 
#'   distribution = pnorm,
#'   log10 = TRUE,
#'   points_params = list(cex = .5)
#' )
#' 
#' # Make same plot with difference between observed and expected values on the y-axis 
#' pp_conf_plot(
#'   obs=smp, 
#'   distribution = pnorm,
#'   difference = TRUE
#' )
#' 
#' # Make same plot with samples plotted as a blue line, expected value line plotted as a red line,
#' # and pointwise bounds plotted as black lines
#' pp_conf_plot(
#'   obs=smp, 
#'   distribution = pnorm,
#'   plot_pointwise = TRUE,
#'   method = "ell",
#'   points_params = list(col="blue", type="l"),
#'   line_params = list(col="red")
#' )
#' 

pp_conf_plot <- function(obs,
                         distribution = pnorm,
                         method = c("ell", "ks"),
                         alpha = 0.05,
                         difference = FALSE,
                         log10 = FALSE,
                         right_tail = FALSE,
                         add = FALSE,
                         dparams = list(),
                         bounds_params = list(),
                         line_params = list(),
                         plot_pointwise = FALSE,
                         pointwise_lines_params = list(),
                         points_params = list(),
                         polygon_params = list(border = NA, col = 'gray'),
                         ...) {
  
  dist_name <- as.character(substitute(distribution))
  
  if(length(dparams) == 0) {
    # equivalence between base R and MASS::fitdistr distribution names
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
    
    # initial value for some distributions
    initVal <- function(distributionName) {
      switch(
        distributionName,
        pbeta = list(shape1 = 1, shape2 = 1),
        pchisq = list(df = 1),
        pf = list(df1 = 1, df2 = 2),
        pt = list(df = 1),
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
            
          } else {
            
            dparams <- MASS::fitdistr(x = obs, densfun = corresp(dist_name))$estimate
            
          }
          
        } else {
          dparams <- MASS::fitdistr(x = obs, densfun = corresp(dist_name), start = initVal(dist_name))$estimate
        }
      }
    })
  }
  
  dots <- list(...)
  method <- match.arg(method)
  if (is.null(dots$ylab)) {
    if (difference && log10 && right_tail) {
      ylab <- expression("-log"[10]*"(1 - Observed probabilities) + log"[10]*"(1 - Expected probabilities)")
    } else if (difference && log10) {
      ylab <- expression("-log"[10]*"(Observed probabilities) + log"[10]*"(Expected probabilities)")
    } else if(difference) {
      ylab <- 'Obseved probabilities - Expected probabilities'
    } else if (log10 && right_tail) {
      ylab <- expression("-log"[10]*"(1 - Observed probabilities)")
    } else if (log10) {
      ylab <- expression("-log"[10]*"(Observed probabilities)")
    } else {
      ylab <- 'Observed probabilities'
    }
  } else {
    ylab <- dots$ylab
    dots <- within(dots, rm(ylab))
  }
  if (is.null(dots$xlab)) {
    if (log10 && right_tail) {
      xlab <- expression("-log"[10]*"(1 - Expected probabilities)")
    } else if (log10) {
      xlab <- expression("-log"[10]*"(Expected probabilities)")
    } else {
      xlab <- 'Expected probabilities'
    }
  } else {
    xlab <- dots$xlab
    dots <- within(dots, rm(xlab))
  }
  
  samp.size <- length(obs)
  conf.int <- 1 - alpha
  conf <- c(alpha / 2, conf.int + alpha / 2)
  ## The observed and expected probabilities. Expected probabilities are based on the specified
  ## distribution
  # constant for visual expansion of confidence regions
  c <- .5 
  obs.pts <- do.call(distribution, c(list(q=sort(obs)), dparams))
  exp.pts <- ppoints(samp.size, a=0)
  if (log10 == TRUE && right_tail == TRUE) {
    
    exp.pts <- -log10(1 - exp.pts)
    low_exp_pt <- c * -log10(1 - (1 / max(samp.size * 1.25, samp.size + 2))) + (1 - c) * exp.pts[1]
    high_exp_pt <- c * -log10((1 / max(samp.size * 1.25, samp.size + 2))) + (1 - c) * exp.pts[samp.size]
    obs.pts <- -log10(1 - obs.pts)
    
  }
  else if (log10 == TRUE) {
    
    exp.pts <- -log10(exp.pts)
    low_exp_pt <- c * -log10(1 / max(samp.size * 1.25, samp.size + 2)) + (1 - c) * exp.pts[1]
    high_exp_pt <- c * -log10(1 - (1 / max(samp.size * 1.25, samp.size + 2))) + (1 - c) * exp.pts[samp.size]
    obs.pts <- -log10(obs.pts)
    
  }
  else {
    
    low_exp_pt <- c * (1 / max(samp.size * 1.25, samp.size + 2)) + (1 - c) * exp.pts[1]
    high_exp_pt <- c * (1 - (1 / max(samp.size * 1.25, samp.size + 2))) + (1 - c) * exp.pts[samp.size]

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
    
    pointwise.low <- qbeta(conf[1], 1:samp.size, samp.size:1)
    pointwise.high <- qbeta(conf[2], 1:samp.size, samp.size:1)
    
    global.bounds <- do.call(get_bounds_two_sided,
                             c(list(alpha = alpha, n = samp.size), bounds_params))

    if (method == "ell") {
      
      global.low <- global.bounds$lower_bound
      global.high <- global.bounds$upper_bound
      
    } else if (method == "ks") {
      
      probs <- ppoints(samp.size)
      epsilon <- sqrt((1 / (2 * samp.size)) * log(2 / (1 - conf.int)))
      global.low <- pmax(probs - epsilon, rep(0, samp.size))
      global.high <- pmin(probs + epsilon, rep(1, samp.size))
      
    }
    
    if (log10 == TRUE && right_tail == TRUE) {
      pointwise.low <- -log10(1 - pointwise.low)
      pointwise.high <- -log10(1 - pointwise.high)
      global.low <- -log10(1 - global.low)
      global.high <- -log10(1 - global.high)
    } else if (log10 == TRUE) {
      pointwise.low <- -log10(pointwise.low)
      pointwise.high <- -log10(pointwise.high)
      global.low <- -log10(global.low)
      global.high <- -log10(global.high)
    }
    
    # code to extend region for visibility
    global.low <- c(global.low[1], global.low, global.low[samp.size])
    global.high <- c(global.high[1], global.high, global.high[samp.size])
    pointwise.low <- c(pointwise.low[1], pointwise.low, pointwise.low[samp.size])
    pointwise.high <- c(pointwise.high[1], pointwise.high, pointwise.high[samp.size])
    exp.pts <- c(low_exp_pt, exp.pts, high_exp_pt)
    
    if ("ylim" %in% names(dots)) {
      
      bottom <- dots$ylim[1] - 1000
      top <- dots$ylim[2] + 1000
      
    } else {
      
      global.low_temp <- global.low[is.finite(global.low)]
      global.high_temp <- global.high[is.finite(global.high)]
      bottom <- min(global.low_temp) - 1000
      top <- max(global.high_temp) + 1000
      
    }
    
    if (difference) {
      do.call(
        polygon, 
        c(list(x = c(exp.pts, rev(exp.pts)),
               y = pmin(pmax(c(global.low - exp.pts, rev(global.high) - rev(exp.pts)), bottom), top)),
          polygon_params)
      )
      if (plot_pointwise) {
        
        do.call(lines, c(list(x = exp.pts, y = pointwise.low - exp.pts), pointwise_lines_params))
        do.call(lines, c(list(x = exp.pts, y = pointwise.high - exp.pts), pointwise_lines_params))
        
      }
      
    } else {
      
      do.call(
        polygon,
        c(list(x = c(exp.pts, rev(exp.pts)),
               y = pmin(pmax(c(global.low, rev(global.high)), bottom), top)), 
          polygon_params)
      )
      if (plot_pointwise) {
        
        do.call(lines, c(list(x = exp.pts, y = pointwise.low), pointwise_lines_params))
        do.call(lines, c(list(x = exp.pts, y = pointwise.high), pointwise_lines_params))
        
      }
      
    }
  }
  do.call(points, c(list(x = exp.pts[2:(samp.size + 1)], y = y.pts), points_params))
  if (difference) {
    
    do.call(lines, c(list(x = c(min(exp.pts), max(exp.pts)), y = c(0, 0)), line_params))
    
  } else{
    
    do.call(lines, c(list(x = c(min(exp.pts), max(exp.pts)), y = c(min(exp.pts), max(exp.pts))), line_params))
    
  }
  
}
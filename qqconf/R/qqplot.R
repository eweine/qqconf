#' QQ Plot with Simultaneous and Pointwise Testing Bounds.
#'
#' Create a qq-plot with with a shaded simultaneous acceptance region and,
#' optionally, lines for a point-wise region. The observed values are plotted
#' against their expected values had they come from the specified distribution.
#'
#' If any of the points of the qq-plot fall outside the simultaneous acceptance region for the selected
#' level alpha test, that means that we can reject the null hypothesis that the data are i.i.d. draws from the
#' specified distribution. If \code{difference} is set to TRUE, the vertical axis plots the
#' observed quantile minus expected quantile. If pointwise bounds are used, then on average, alpha * n of the points will fall outside
#' the bounds under the null hypothesis, so the chance that the qq-plot has any points falling outside of the pointwise bounds
#' is typically much higher than alpha under the null hypothesis. For this reason, a simultaneous region is preferred.
#'
#' @param obs The observed data.
#' @param distribution The quantile function for the specified distribution. Defaults to \code{qnorm}.
#' Custom distributions are allowed as long as all parameters are supplied in dparams.
#' @param method Method for simultaneous testing bands. Must be either "ell" (equal local levels test), which applies a level \eqn{\eta} pointwise
#' test to each order statistic such that the Type I error of the global test is \code{alpha}, or "ks" to apply a
#' Kolmogorov-Smirnov test. "ell" is recommended.
#' @param alpha Type I error of global test of whether the data come from the reference distribution.
#' @param difference Whether to plot the difference between the observed and
#'   expected values on the vertical axis.
#' @param log10 Whether to plot axes on -log10 scale (e.g. to see small p-values). Can only be used for strictly
#' positive distributions.
#' @param right_tail This argument is only used if \code{log10} is \code{TRUE}. When \code{TRUE},
#' the x-axis is -log10(1 - Expected Quantile) and the y-axis is -log10(1 - Observed Quantile).
#' When \code{FALSE} (default) the x-axis is -log10(Expected Quantile) and the y-axis is
#' -log10(Observed Quantile). The argument should be set to \code{TRUE} only when the support
#' of the distribution lies in (0, 1), and one wants to make
#' observations in the right tail of the distribution easier to see. The argument should be
#' set to \code{FALSE} when one wants to make
#' observations in the left tail of the distribution easier to see.
#' @param add Whether to add points to an existing plot.
#' @param dparams List of additional arguments for the quantile function of the distribution
#'   (e.g. df=1). Note that if any parameters of the distribution are specified, parameter estimation will not be performed
#'   on the unspecified parameters, and instead they will take on the default values set by the distribution function.
#'   For the uniform distribution, parameter estimation is not performed, and
#'   the default parameters are max = 1 and min = 0.
#'   For other distributions parameters will be estimated if not provided.
#'   For the normal distribution, we estimate the mean as the median and the standard deviation as \eqn{Sn} from the paper by Rousseeuw and Croux 1993
#'   "Alternatives to the Median Absolute Deviation". For all other distributions besides uniform and normal,
#'   the code uses MLE to estimate the parameters. Note that estimation is not implemented for custom distributions, so all
#'   parameters of the distribution must be provided by the user.
#' @param bounds_params List of optional arguments for \code{get_bounds_two_sided}
#'   (i.e. \code{tol}, \code{max_it}, \code{method}).
#' @param line_params Arguments passed to the \code{lines} function to modify the line that indicates a perfect fit of the
#'   reference distribution.
#' @param plot_pointwise Boolean indicating whether pointwise bounds should be added to the plot
#' @param pointwise_lines_params Arguments passed to the \code{lines} function that modifies pointwise bounds when plot_pointwise is
#'   set to TRUE.
#' @param points_params Arguments to be passed to the \code{points} function to plot the data.
#' @param polygon_params Arguments to be passed to the polygon function to construct simultaneous confidence region.
#'   By default \code{border} is set to NA and \code{col} is set to grey.
#' @param prob_pts_method (optional) method used to get probability points for plotting.
#' The quantile function will be applied to these points to
#' get the expected values.  When this argument is set to \code{"normal"}
#' (recommended for a normal QQ plot) \code{ppoints(n)} will be used,  which is what
#' most other plotting software uses. When this argument is set to \code{"uniform"}
#' (recommended for a uniform QQ plot) \code{ppoints(n, a=0)}, which are the expected
#' values of the order statistics of Uniform(0, 1), will be used.  Finally,
#'  when this argument is set to \code{"median"} (recommended for all other
#'  distributions) \code{qbeta(.5, c(1:n), c(n:1))} will be used. Under the default
#'  setting, \code{"best_available"}, the probability points as recommended above will
#'  be used. Note that \code{"median"} is suitable for all distributions and is
#'  particularly recommended when alpha is large.
#' @param ... Additional arguments passed to the plot function.
#'
#' @importFrom stats median pnorm ppoints qbeta qnorm
#' @importFrom graphics abline lines plot points polygon
#'
#' @export
#'
#' @return None, QQ plot is produced.
#'
#' @examples
#' set.seed(0)
#' smp <- runif(100)
#'
#' # Plot QQ plot against uniform(0, 1) distribution
#' qq_conf_plot(
#'   obs=smp,
#'   distribution = qunif
#' )
#'
#' # Make same plot on -log10 scale to highlight small p-values,
#' # with radius of plot circles also reduced by .5
#' qq_conf_plot(
#'   obs=smp,
#'   distribution = qunif,
#'   points_params = list(cex = .5),
#'   log10 = TRUE
#' )
#'
#' # Make same plot with difference between observed and expected values on the y-axis
#' qq_conf_plot(
#'   obs=smp,
#'   distribution = qunif,
#'   difference = TRUE
#' )
#'
#' # Make same plot with sample plotted as a blue line, expected value line plotted as a red line,
#' # and with pointwise bounds plotted as black lines
#' qq_conf_plot(
#'   obs=smp,
#'   distribution = qunif,
#'   plot_pointwise = TRUE,
#'   points_params = list(col="blue", type="l"),
#'   line_params = list(col="red")
#' )
#'

qq_conf_plot <- function(obs,
                         distribution = qnorm,
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
                         prob_pts_method = c("best_available", "normal", "uniform", "median"),
                         ...) {

  if (!("p" %in% names(formals(distribution)))) {

    stop("distribution function must take 'p' as an argument.
         Did you mean to make a PP plot?")

  }

  dots <- list(...)
  method <- match.arg(method)
  if ( is.null(dots$ylab)) {
    if (difference && log10 && right_tail) {
      ylab <- expression("-log"[10]*"(1 - Observed quantiles) + log"[10]*"(1 - Expected quantiles)")
    } else if(difference && log10) {
      ylab <- expression("-log"[10]*"(Observed quantiles) + log"[10]*"(Expected quantiles)")
    } else if(difference) {
      ylab <- 'Obseved quantiles - Expected quantiles'
    } else if(log10 && right_tail) {
      ylab <- expression("-log"[10]*"(1 - Observed quantiles)")
    }else if (log10) {
      ylab <- expression("-log"[10]*"(Observed quantiles)")
    } else {
      ylab <- 'Observed quantiles'
    }
  } else {
    ylab <- dots$ylab
    dots <- within(dots, rm(ylab))
  }
  if (is.null(dots$xlab)) {
    if (log10 && right_tail) {
      xlab <- expression("-log"[10]*"(1 - Expected quantiles)")
    } else if(log10) {
      xlab <- expression("-log"[10]*"(Expected quantiles)")
    }
    else {
      xlab <- 'Expected quantiles'
    }
  } else {
    xlab <- dots$xlab
    dots <- within(dots, rm(xlab))
  }

  samp_size <- length(obs)
  conf_int <- 1 - alpha
  conf <- c(alpha / 2, conf_int + alpha / 2)
  ## The observed and expected quantiles. Expected quantiles are based on the specified
  ## distribution
  # constant for visual expansion of confidence regions
  c <- .5
  obs_pts <- sort(obs)

  dist_name <- as.character(substitute(distribution))
  if(length(dparams) == 0 && dist_name == "qunif") {

    dparams['min'] <- 0
    dparams['max'] <- 1

  } else if (length(dparams) == 0) {

    cat("no dparams supplied. Estimating parameters from the data...\n")
    MASS_name <- get_mass_name_from_distr(dist_name, "qq")
    dparams <- estimate_params_from_data(MASS_name, obs)

  }

  prob_pts_method <- match.arg(prob_pts_method)

  if (prob_pts_method == "best_available") {

    prob_pts_method <- get_best_available_prob_pts_method(dist_name)

  }

  global_bounds <- get_qq_band(
    obs = obs,
    alpha = alpha,
    distribution = distribution,
    dparams = dparams,
    ell_params = bounds_params,
    band_method = method,
    prob_pts_method = prob_pts_method
  )

  global_low <- global_bounds$lower_bound
  global_high <- global_bounds$upper_bound

  exp_pts <- global_bounds$expected_value

  ext_quantile <- get_extended_quantile(prob_pts_method, samp_size)

  if (log10 && right_tail) {

    exp_pts <- -log10(1 - exp_pts)
    low_exp_pt <- c * -log10(do.call(distribution, c(list(p=ext_quantile$high_pt), dparams))) + (1 - c) * exp_pts[1]
    high_exp_pt <- c * -log10(do.call(distribution, c(list(p=ext_quantile$low_pt), dparams))) + (1 - c) * exp_pts[samp_size]
    if (any(obs_pts <= 0)) {

      stop("log10 scaling can only be used with strictly positive distributions. Consider using pp_conf_plot.")

    }
    obs_pts <- -log10(1 - obs_pts)

  }
  else if (log10 == TRUE) {

    exp_pts <- -log10(exp_pts)
    low_exp_pt <- c * -log10(do.call(distribution, c(list(p=ext_quantile$low_pt), dparams))) + (1 - c) * exp_pts[1]
    high_exp_pt <- c * -log10(do.call(distribution, c(list(p=ext_quantile$high_pt), dparams))) + (1 - c) * exp_pts[samp_size]
    if (any(obs_pts <= 0)) {

      stop("log10 scaling can only be used with strictly positive distributions. Consider using pp_conf_plot.")

    }
    obs_pts <- -log10(obs_pts)

  }
  else {

    low_exp_pt <- c * do.call(distribution, c(list(p=ext_quantile$low_pt), dparams)) + (1 - c) * exp_pts[1]
    high_exp_pt <- c * do.call(distribution, c(list(p=ext_quantile$high_pt), dparams)) + (1 - c) * exp_pts[samp_size]

  }
  if (difference) {
    y_pts <- obs_pts - exp_pts
  } else {
    y_pts <- obs_pts
  }

  ## When not adding points to a qq-plot compute pointwise and global confidence bounds.
  if (!add) {
    left <- exp_pts[1]
    right <- exp_pts[samp_size]
    bottom <- min(y_pts) #obs_pts[1]
    top <- max(y_pts) #obs_pts[samp_size]
    do.call(plot, c(list(x=c(left, right), y=c(bottom, top), type='n', xlab=xlab, ylab=ylab), dots))
    pointwise_low <- do.call(distribution,
                             c(list(p=qbeta(conf[1], 1:samp_size, samp_size:1)), dparams))
    pointwise_high <- do.call(distribution,
                              c(list(p=qbeta(conf[2], 1:samp_size, samp_size:1)), dparams))

    if(log10 == TRUE && right_tail == TRUE) {

      pointwise_low <- -log10(1 - pointwise_low)
      pointwise_high <- -log10(1 - pointwise_high)
      global_low <- -log10(1 - global_low)
      global_high <- -log10(1 - global_high)

    }
    else if (log10 == TRUE) {

      pointwise_low <- -log10(pointwise_low)
      pointwise_high <- -log10(pointwise_high)
      global_low <- -log10(global_low)
      global_high <- -log10(global_high)

    }

    if ("ylim" %in% names(dots)) {

      bottom <- dots$ylim[1] - 1000
      top <- dots$ylim[2] + 1000

    } else {

      global_low_temp <- global_low[is.finite(global_low)]
      global_high_temp <- global_high[is.finite(global_high)]
      bottom <- min(global_low_temp) - 1000
      top <- max(global_high_temp) + 1000

    }

    if (difference) {

      low_global_diff <- global_low - exp_pts
      low_global_diff <- c(low_global_diff[1], low_global_diff, low_global_diff[samp_size])
      high_global_diff <- global_high - exp_pts
      high_global_diff <- c(high_global_diff[1], high_global_diff, high_global_diff[samp_size])
      low_pointwise_diff <- pointwise_low - exp_pts
      low_pointwise_diff <- c(low_pointwise_diff[1], low_pointwise_diff, low_pointwise_diff[samp_size])
      high_pointwise_diff <- pointwise_high - exp_pts
      high_pointwise_diff <- c(high_pointwise_diff[1], high_pointwise_diff, high_pointwise_diff[samp_size])
      exp_pts <- c(low_exp_pt, exp_pts, high_exp_pt)

      do.call(
        polygon,
        c(list(x = c(exp_pts, rev(exp_pts)),
               y = pmin(pmax(c(low_global_diff, rev(high_global_diff)), bottom), top)),
          polygon_params)
      )
      if (plot_pointwise) {

        do.call(lines, c(list(x = exp_pts, y = low_pointwise_diff), pointwise_lines_params))
        do.call(lines, c(list(x = exp_pts, y = high_pointwise_diff), pointwise_lines_params))

      }

    } else {

      # code to extend region for visibility
      global_low <- c(global_low[1], global_low, global_low[samp_size])
      global_high <- c(global_high[1], global_high, global_high[samp_size])
      pointwise_low <- c(pointwise_low[1], pointwise_low, pointwise_low[samp_size])
      pointwise_high <- c(pointwise_high[1], pointwise_high, pointwise_high[samp_size])
      exp_pts <- c(low_exp_pt, exp_pts, high_exp_pt)

      do.call(
        polygon,
        c(list(x = c(exp_pts, rev(exp_pts)),
               y = pmin(pmax(c(global_low, rev(global_high)), bottom), top)),
          polygon_params)
        )
      if (plot_pointwise) {

        do.call(lines, c(list(x = exp_pts, y = pointwise_low), pointwise_lines_params))
        do.call(lines, c(list(x = exp_pts, y = pointwise_high), pointwise_lines_params))

      }

    }

    do.call(points, c(list(x = exp_pts[2:(samp_size + 1)], y = y_pts), points_params))

  } else {

    do.call(points, c(list(x = exp_pts, y = y_pts), points_params))

  }

  if (difference) {

    do.call(lines, c(list(x = c(min(exp_pts), max(exp_pts)), y = c(0, 0)), line_params))

  } else{

    do.call(lines, c(list(x = c(min(exp_pts), max(exp_pts)), y = c(min(exp_pts), max(exp_pts))), line_params))

  }

}

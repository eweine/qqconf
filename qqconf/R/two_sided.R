#' Check Validity of Two-Sided Bounds
#'
#' Given bounds for a two sided test, this checks that none of
#' the bounds fall outside of [0, 1] and that all upper bounds
#' are greater than the corresponding lower bounds.
#' This also ensures the the length of the bounds are the same.
#' This not meant to be called by the user.
#'
#' @param lower_bounds Numeric vector where the ith component is the lower bound
#' for the ith order statistic.
#' @param upper_bounds Numeric vector where the ith component is the lower bound
#' for the ith order statistic.
#'
#' @return None
#'
check_bounds_two_sided <- function(lower_bounds,
                                   upper_bounds) {

  if(any(lower_bounds > 1) || any(lower_bounds < 0)) {

    stop("Not all lower bounds between 0 and 1 (inclusive)")

  }

  if(any(upper_bounds > 1) || any(upper_bounds < 0)) {

    stop("Not all upper bounds between 0 and 1 (inclusive)")

  }

  if(any(upper_bounds - lower_bounds < 0)) {

    stop("Not all upper bounds are greater than their corresponding lower bounds")

  }

  if(length(lower_bounds) != length(upper_bounds)) {

    stop("The length of the bounds differ")

  }

  if(is.unsorted(lower_bounds)) {

    stop("Only lower bounds in ascending order are supported")

  }

  if(is.unsorted(upper_bounds)) {

    stop("Only upper bounds in ascending order are supported")

  }

}


#' Monte Carlo Simulation for Two-Sided Test
#'
#' Given bounds for a two sided test on uniform order statistics, this computes
#' the Type I Error Rate \eqn{\alpha} using simulations.
#'
#' @param lower_bounds Numeric vector where the ith component is the lower bound
#' for the ith order statistic. The components must be distinct values in (0, 1) that
#' are in ascending order.
#' @param upper_bounds Numeric vector where the ith component is the lower bound
#' for the ith order statistic. The values must be in ascending order and
#' the ith component must be larger than the ith component of the lower bounds.
#' @param num_sims (optional) Number of simulations to be run, 1 Million by default.
#'
#' @return Type I Error Rate \eqn{\alpha}
#'
monte_carlo_two_sided <- function(lower_bounds,
                                  upper_bounds,
                                  num_sims = 1000000) {

  check_bounds_two_sided(lower_bounds, upper_bounds)
  n <- length(lower_bounds)

  num_in_ci <- 0 # tracks how many sims have all points in the confidence intervals

  for(iter in seq(from = 1, to = num_sims, by = 1)) {

    # generate uniform RVs
    iter_rvs <- stats::runif(n)
    # sort them
    iter_order_stats <- sort(iter_rvs)

    if (all(iter_order_stats > lower_bounds) && all(iter_order_stats < upper_bounds)) {

      num_in_ci <- num_in_ci + 1

    }

  }

  alpha <- 1 - (num_in_ci / num_sims)
  return(alpha)

}

#' Calculates Approximate Local Level
#'
#' This function uses the approximation from Gontscharuk & Finner's Asymptotics of
#' goodness-of-fit tests based on minimum p-value statistics (2017) to approximate
#' local levels for finite sample size. We use these authors constants for \eqn{\alpha} = .1, and .05,
#' and for \eqn{\alpha} = .01 we use a slightly different approximation.
#'
#' @param n Number of tests to do
#' @param alpha Global type I error rate \eqn{\alpha} of the tests
#'
#' @return Approximate local level
get_asymptotic_approx_corrected_alpha <- function(n, alpha) {

  if (alpha == .01) {

    c_alpha <- 1.591

  } else if (alpha == .05) {

    c_alpha <- 1.3

  } else if (alpha == .1) {

    c_alpha <- 1.1

  }

  eta_approx = -log(1-alpha)/(2*log(log(n))*log(n))*(1-c_alpha*log(log(log(n)))/log(log(n)))
  return(eta_approx)

}

#' Calculates Global Significance Level From Simultaneous Two-Sided Bounds for Rejection Region
#'
#' For a test of uniformity of i.i.d. observations on the unit interval, this function will determine the significance
#' level as a function of the rejection region. Suppose \eqn{n} observations are drawn i.i.d. from some CDF F(x) on the unit interval,
#' and it is desired to test the null hypothesis that F(x) = x for all x in (0, 1) against a two-sided alternative.
#' Suppose the acceptance region for the test is described by a set of intervals, one for each order statistic.
#' Given the bounds for these intervals, this function calculates the significance level of the test where the
#' null hypothesis is rejected if at least one of the order statistics is outside its corresponding interval.
#'
#' Uses the method of Moscovich and Nadler (2016) as implemented in Crossprob (Moscovich 2020).
#'
#' @param lower_bounds Numeric vector where the ith component is the lower bound for the acceptance interval
#' for the ith order statistic. The components must lie in [0, 1], and each component must be greater than
#' or equal to the previous one.
#' @param upper_bounds Numeric vector of the same length as \code{lower_bounds} where the ith component is the upper bound
#' for the acceptance interval for the ith order statistic. The components must lie in [0, 1], and each component must be
#' greater than or equal to the previous one. In addition,
#' the ith component of \code{upper_bounds} must be greater than or equal to the ith component of \code{lower_bounds}.
#'
#' @return Global Significance Level \eqn{\alpha}
#'
#' @examples
#' # For X1, X2 iid unif(0,1), calculate 1 - P(.1 < min(X1, X2) < .6 and .5 < max(X1, X2) < .9)
#' get_level_from_bounds_two_sided(lower_bounds = c(.1, .5), upper_bounds = c(.6, .9))
#'
#' # Finds the global significance level corresponding to the local level eta.
#' # Suppose we reject the null hypothesis that X1, ..., Xn are iid unif(0, 1) if and only if at least
#' # one of the order statistics X(i) is significantly different from
#' # its null distribution based on a level-eta
#' # two-sided test, i.e. we reject if and only if X(i) is outside the interval
#' # (qbeta(eta/2, i, n - i + 1), qbeta(1 - eta/2, i, n - i + 1)) for at least one i.
#' # The lines of code below calculate the global significance level of
#' # the test (which is necessarily larger than eta if n > 1).
#' n <- 100
#' eta <- .05
#' lb <- qbeta(eta / 2, c(1:n), c(n:1))
#' ub <- qbeta(1 - eta / 2, c(1:n), c(n:1))
#' get_level_from_bounds_two_sided(lower_bounds = lb, upper_bounds = ub)
#'
#' @references
#' \itemize{
#' \item{\href{https://www.sciencedirect.com/science/article/abs/pii/S0167715216302802}{
#' Moscovich, Amit, and Boaz Nadler. "Fast calculation of boundary crossing probabilities for Poisson processes."
#' Statistics & Probability Letters 123 (2017): 177-182.}}
#' \item{\href{https://github.com/mosco/crossing-probability}{
#' Amit Moscovich (2020). Fast calculation of p-values for one-sided
#' Kolmogorov-Smirnov type statistics. arXiv:2009.04954}}
#' }
#'
#' @useDynLib qqconf
#'
#' @export
get_level_from_bounds_two_sided <- function(lower_bounds,
                                            upper_bounds) {

  check_bounds_two_sided(lower_bounds, upper_bounds)
  alpha <- 1 - fft_get_level_from_bounds_two_sided(lower_bounds, upper_bounds)
  return(alpha)

}

#' Calculates Rejection Region of Two-Sided Equal Local Levels Test.
#'
#' The context is that n i.i.d. observations are assumed to be drawn
#' from some distribution on the unit interval with c.d.f. F(x), and it is
#' desired to test the null hypothesis that F(x) = x for all x in (0,1),
#' referred to as the "global null hypothesis," against a two-sided alternative.
#' An "equal local levels" test is used, in which each of the n order statistics is
#' tested for significant deviation from its null distribution by a 2-sided test
#' with significance level \eqn{\eta}.  The global null hypothesis is rejected if at
#' least one of the order statistic tests is rejected at level \eqn{\eta}, where \eqn{\eta} is
#' chosen so that the significance level of the global test is alpha.
#' Given the size of the dataset n and the desired global significance level alpha,
#' this function calculates the local level \eqn{\eta} and the acceptance/rejection regions for the test.
#' There are a set of n intervals, one for each order statistic.
#' If at least one order statistic falls outside the corresponding interval,
#' the global test is rejected.
#'
#'
#' @param alpha Desired global significance level of the test.
#' @param n Size of the dataset.
#' @param tol (optional) Relative tolerance of the \code{alpha} level of the
#' simultaneous test. Defaults to 1e-8. Used only if \code{method} is set to
#' "search" or if method is set to "best_available" and the best available
#' method is a search.
#' @param max_it (optional) Maximum number of iterations of Binary Search Algorithm
#' used to find the bounds. Defaults to 100 which should be much larger than necessary
#' for a reasonable tolerance. Used only if \code{method} is set to
#' "search" or if method is set to "best_available" and the best available
#' method is a search.
#' @param method (optional) Argument indicating if the calculation should be done using a highly
#' accurate approximation, "approximate", or if the calculations should be done using an exact
#' binary search calculation, "search". The default is "best_available" (recommended), which uses the exact search
#' when either (i) the approximation isn't available or (ii) the approximation is available but isn't highly accurate and the search method
#' isn't prohibitively slow (occurs for small to moderate \code{n} with \code{alpha} = .1).
#' Of note, the approximate method is only available for alpha values of .1, .05, and .01. In the case of alpha = .05 or .01, the
#' approximation is highly accurate for all values of \code{n} up to at least \code{10^6}.
#'
#' @return A list with components
#' \itemize{
#'   \item lower_bound - Numeric vector of length \code{n} containing the lower bounds
#'   for the acceptance regions of the test of each order statistic.
#'   \item upper_bound - Numeric vector of length \code{n} containing the upper bounds
#'   for the acceptance regions of the test of each order statistic.
#'   \item x - Numeric vector of length \code{n} containing the expectation of each order statistic.
#'   These are the x-coordinates for the bounds if used in a pp-plot.
#'   The value is \code{c(1:n) / (n + 1)}.
#'   \item local_level - Significance level \eqn{\eta} of the local test on each individual order statistic. It is equal for all order
#'   statistics and will be less than \code{alpha} for all \code{n} > 1.
#' }
#'
#' @examples
#' get_bounds_two_sided(alpha = .05, n = 100)
#'
#' @importFrom utils head tail
#'
#' @export
get_bounds_two_sided <- function(alpha,
                                n,
                                tol = 1e-8,
                                max_it = 100,
                                method=c("best_available", "approximate", "search")) {

  if (alpha >= 1 || alpha <= 0) {

    stop("alpha must be between 0 and 1 (exclusive).")

  }

  if (as.integer(n) != n) {

    stop("n must be an integer")

  }

  if (n < 1) {

    stop("n must be greater than 0")

  }

  if (n == 1) {

    return(list(lower_bound = c(alpha / 2),
                upper_bound = c(1 - alpha / 2),
                x = c(.5),
                local_level = alpha))

  }

  n_param <- n
  if (n >= 10) {

    # Approximation only available for n < 10
    method <- match.arg(method)

  } else {

    method <- "search"

  }

  # Value used for testing if alpha is approximately equal to a set of pre-set values
  alpha_epsilon <- 10 ^ (-5)

  # Approximations are only available for alpha = .05 or alpha = .01
  if (method == "search") {

    eta_high <- alpha
    eta_low <- alpha / n
    eta_curr <- eta_low + (eta_high - eta_low) / 2
    n_it <- 0

    while (n_it < max_it) {

      n_it <- n_it + 1
      h_vals <- stats::qbeta(eta_curr / 2, 1:n, n:1)
      g_vals <- stats::qbeta(1 - (eta_curr / 2), 1:n, n:1)
      test_alpha <- get_level_from_bounds_two_sided(h_vals, g_vals)

      if (abs(test_alpha - alpha) / alpha <= tol) break

      if (test_alpha > alpha) {

        eta_high <- eta_curr
        eta_curr <- eta_curr - (eta_curr - eta_low) / 2

      } else if (test_alpha < alpha) {

        eta_low <- eta_curr
        eta_curr <- eta_curr + (eta_high - eta_curr) / 2

      }

      eta <- eta_curr

    }

    if(n_it == max_it) {

      warning("Maximum number of iterations reached.")

    }

  }

  else if (method == "approximate") {

    if (!(between(alpha, .01 - alpha_epsilon, .01 + alpha_epsilon) ||
          between(alpha, .05 - alpha_epsilon, .05 + alpha_epsilon) ||
          between(alpha, .1 - alpha_epsilon, .1 + alpha_epsilon))) {

      stop("The approximate method is only configured for alpha = .1, .05, .01. Consider setting method='search'")

    }

    if (between(alpha, .01 - alpha_epsilon, .01 + alpha_epsilon)) {

      lookup_table <- alpha_01_df

    } else if (between(alpha, .05 - alpha_epsilon, .05 + alpha_epsilon)) {

      lookup_table <- alpha_05_df

    }

    if (
      between(alpha, .01 - alpha_epsilon, .01 + alpha_epsilon) ||
        between(alpha, .05 - alpha_epsilon, .05 + alpha_epsilon)
        ) {

      if (n %in% lookup_table$n) {

        eta_df <- lookup_table[lookup_table$n == n_param, ]
        eta <- eta_df$local_level[1]

      } else if (
        (between(alpha, .05 - alpha_epsilon, .05 + alpha_epsilon) && n > 10 ^ 5) ||
          (between(alpha, .01 - alpha_epsilon, .01 + alpha_epsilon) && n > 10 ^ 5)) {

        eta <- get_asymptotic_approx_corrected_alpha(n, alpha)

      }
      else {

        # Do linear interpolation
        larger_n_df <- lookup_table[lookup_table$n > n_param, ]
        larger_n_df <- head(larger_n_df, 1)
        larger_n <- larger_n_df$n[1]
        larger_n_eta <- larger_n_df$local_level[1]

        smaller_n_df <- lookup_table[lookup_table$n < n_param, ]
        smaller_n_df <- tail(smaller_n_df, 1)
        smaller_n <- smaller_n_df$n[1]
        smaller_n_eta <- smaller_n_df$local_level[1]

        # y = mx + b
        m <- (larger_n_eta - smaller_n_eta) / (larger_n - smaller_n)
        b <- smaller_n_eta - m * smaller_n
        eta <- m * n + b

      }

    } else {

      eta <- get_asymptotic_approx_corrected_alpha(n, alpha)

    }

    h_vals <- stats::qbeta(eta / 2, 1:n, n:1)
    g_vals <- stats::qbeta(1 - (eta / 2), 1:n, n:1)

  }

  else if (method == "best_available") {

    # Don't return approximation for alpha = .1 and n <= 500, it's not accurate
    if (between(alpha, .1 - alpha_epsilon, .1 + alpha_epsilon) && n <= 500 ||
        !(between(alpha, .01 - alpha_epsilon, .01 + alpha_epsilon) ||
          between(alpha, .05 - alpha_epsilon, .05 + alpha_epsilon) ||
          between(alpha, .1 - alpha_epsilon, .1 + alpha_epsilon))
        ) {

      return(
        get_bounds_two_sided(
          alpha = alpha,
          n = n,
          tol = tol,
          max_it = max_it,
          method="search"
        )
      )

    } else {

      return(
        get_bounds_two_sided(
          alpha = alpha,
          n = n,
          tol = tol,
          max_it = max_it,
          method="approximate"
        )
      )

    }

  }

  alpha_vec <- seq(from = 1, to = n, by = 1)
  beta_vec <- n - alpha_vec + 1
  order_stats_mean <- alpha_vec / (alpha_vec + beta_vec)

  return(list(lower_bound = h_vals,
              upper_bound = g_vals,
              x = order_stats_mean,
              local_level = eta))

}

#' Create QQ Plot Testing Band
#'
#' Flexible interface for creating a testing band for a Quantile-Quantile (QQ)
#' plot.
#'
#' @param n,obs either a number of observations (specified by setting \code{n}),
#' or a numeric vector of observations (specified by setting \code{obs}). One
#' argument must be specified. If all parameters of \code{distribution}
#' are known, then the testing band only depends on the number of observations \code{n}.
#' Thus, providing \code{n} is simpler when all parameters
#' of \code{distribution} are known and specified via \code{dparams}
#' (or when using all default parameter choices of \code{distribution} is desired).
#' If estimating parameters from the data is preferred, \code{obs} should be
#' specified and estimation will take place as described in the documentation for
#' argument \code{dparams}.
#' @param alpha (optional) desired significance level of the testing band. If \code{method}
#' is set to \code{"ell"} or \code{"ks"}, then this is the global significance
#' level of the testing band. If \code{method} is set to \code{"pointwise"}, then
#' the band is equivalent to simply conducting a level \code{alpha} test on each
#' order statistic individually. Pointwise bands will generally have much
#' larger global Type I error than \code{alpha}. Defaults to \code{.05}.
#' @param distribution The quantile function for the specified distribution.
#' Defaults to \code{qnorm}, which is appropriate for testing normality
#' of the observations in a QQ plot.
#' @param dparams (optional)  List of additional arguments for the \code{distribution} function
#'   (e.g. df=1). If \code{obs} is not specified and this argument is left blank,
#'   the default arguments of \code{distribution} are used. If \code{obs} is specified and this argument is left blank,
#'   parameters are estimated from the data (except if \code{distribution} is set to \code{qunif},
#'   in which case no estimation occurs and the default parameters are \code{max = 1} and \code{min = 0}).
#'   For the normal distribution, we estimate the mean as the median and the standard deviation as \eqn{Sn} from the paper by Rousseeuw and Croux 1993
#'   "Alternatives to the Median Absolute Deviation". For all other distributions besides uniform and normal,
#'   the code uses MLE to estimate the parameters. Note that if any parameters of the distribution are specified
#'   in \code{dparams}, parameter estimation will not be performed
#'   on the unspecified parameters, and instead they will take on the default values
#'   set by \code{distribution}.
#' @param ell_params (optional) list of optional arguments for \code{get_bounds_two_sided}
#'   (i.e. \code{tol}, \code{max_it}, \code{method}). Only used if \code{method}
#'   is set to \code{"ell"}
#' @param band_method (optional) method for creating the testing band. The default,
#' \code{"ell"} uses the equal local levels method
#' (see \code{get_bounds_two_sided} for more information). \code{"ks"} uses
#' the Kolmogorov-Smirnov test. \code{"pointwise"} uses a pointwise band (see
#' documentation for argument \code{alpha} for more information). \code{"ell"}
#' is recommended and is the default.
#' @param prob_pts_method (optional) method used to get probability points for
#' use in a QQ plot. The quantile function will be applied to these points to
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
#'
#' @return A list with components
#' \itemize{
#'   \item lower_bound - Numeric vector of length \code{n} containing the lower bounds
#'   for the acceptance regions of the test corresponding to each order statistic.
#'   These form the lower boundary of the testing band for the QQ-plot.
#'   \item upper_bound - Numeric vector of length \code{n} containing the upper bounds
#'   for the acceptance regions of the test corresponding to each order statistic.
#'   These form the upper boundary of the testing band for the QQ-plot.
#'   \item expected_value - Numeric vector of length \code{n} containing the
#'   exact or approximate expectation (or median) of each order statistic, depending on how
#'   \code{prob_pts_method} is set.
#'   These are the x-coordinates for both the bounds and the data points
#'   if used in a qq-plot. Note that
#'   if adding a band to an already existing plot, it is essential that the same
#'   x-coordinates be used for the bounds as were used to plot the data. Thus,
#'   if some other x-coordinates have been used to plot the data those same
#'   x-coordinates should always be used instead of this vector to plot the bounds.
#'   \item dparams - List of arguments used to apply \code{distribution} to
#'   \code{obs} (if observations are provided). If the user provides parameters,
#'   then these parameters will simply be returned. If parameters are estimated
#'   from the data, then the estimated parameters will be returned.
#' }
#' @export
#'
#' @examples
#'
#' # Get ell level .05 QQ testing band for normal(0, 1) distribution with 100 observations
#' band <- get_qq_band(n = 100)
#'
#' # Get ell level .05 QQ testing band for normal distribution with unknown parameters
#' obs <- rnorm(100)
#' band <- get_qq_band(obs = obs)
#'
#' # Get ell level .05 QQ testing band for t(2) distribution with 100 observations
#' band <- get_qq_band(
#'   n = 100, distribution = qt, dparams = list(df = 2)
#' )
#'
get_qq_band <- function(
  n,
  obs,
  alpha = .05,
  distribution = qnorm,
  dparams = list(),
  ell_params = list(),
  band_method = c("ell", "ks", "pointwise"),
  prob_pts_method = c("best_available", "normal", "uniform", "median")
) {

  if (missing(obs) && missing(n)) {

    stop("one of obs or n must be supplied")

  }

  if (!("p" %in% names(formals(distribution)))) {

    stop("distribution function must take 'p' as an argument.")

  }

  band_method <- match.arg(band_method)
  prob_pts_method <- match.arg(prob_pts_method)
  dist_name <- as.character(substitute(distribution))

  if(!missing(obs)) {

    n <- length(obs)
    if (length(dparams) == 0) {

      if (dist_name %in% c("qunif", "punif")) {

        dparams['min'] <- 0
        dparams['max'] <- 1

      }
      else {

        cat("no dparams supplied. Estimating parameters from the data...\n")
        MASS_name <- get_mass_name_from_distr(dist_name, "qq")
        dparams <- estimate_params_from_data(MASS_name, obs)

      }

    }

  }

  if (prob_pts_method == "best_available") {

    prob_pts_method <- get_best_available_prob_pts_method(dist_name)

  }

  if (prob_pts_method == "uniform") {

    raw_exp_pts <- ppoints(n, a=0)

  } else if (prob_pts_method == "normal") {

    raw_exp_pts <- ppoints(n)

  } else if (prob_pts_method == "median") {

    raw_exp_pts <- qbeta(.5, c(1:n), c(n:1))

  }

  tryCatch(
    expr = {

      exp_pts <- do.call(distribution, c(list(p=raw_exp_pts), dparams))

    },
    error = function(e){
      message('Error applying distribution to calculated bounds.')
      message('Did you specify all necessary parameters of distribution?')
      message(e)
    }
  )

  if (band_method == "ell") {

    ell_params["n"] <- n
    ell_params["alpha"] <- alpha
    ell_bounds <- do.call(get_bounds_two_sided, ell_params)
    lower_bound <- ell_bounds$lower_bound
    upper_bound <- ell_bounds$upper_bound

  } else if (band_method == "ks") {

    probs <- ppoints(n)
    epsilon <- sqrt((1 / (2 * n)) * log(2 / alpha))
    lower_bound <- pmax(probs - epsilon, rep(0, n))
    upper_bound <- pmin(probs + epsilon, rep(1, n))

  } else if (band_method == "pointwise") {

    conf.int <- 1 - alpha
    conf <- c(alpha / 2, conf.int + alpha / 2)
    lower_bound <- qbeta(conf[1], 1:n, n:1)
    upper_bound <- qbeta(conf[2], 1:n, n:1)

  }

  lower_bound <- do.call(distribution, c(list(p = lower_bound), dparams))
  upper_bound <- do.call(distribution, c(list(p = upper_bound), dparams))

  return(
    list(
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      expected_value = exp_pts,
      dparams = dparams
    )
  )

}

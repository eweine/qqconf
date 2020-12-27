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

    stop("Not all lower bounds between 0 and 1")

  }

  if(any(upper_bounds > 1) || any(upper_bounds < 0)) {

    stop("Not all upper bounds between 0 and 1")

  }

  if(any(upper_bounds - lower_bounds < 0)) {

    stop("Not all upper bounds are greater than their corresponding lower bounds")

  }

  if(length(lower_bounds) != length(upper_bounds)) {

    stop("The length of the bounds differ")

  }

}


#' Monte Carlo Simulation for Two-Sided Test
#'
#' Given bounds for a two sided test on uniform order statistics, this computes
#' the Type I Error Rate \eqn{\alpha} using simulations.
#'
#' @param lower_bounds Numeric vector where the ith component is the lower bound
#' for the ith order statistic. The values must be in ascending order.
#' @param upper_bounds Numeric vector where the ith component is the lower bound
#' for the ith order statistic. The values must be in ascending order and
#' the ith component must be larger than the ith component of the lower bounds.
#' @param num_sims (Optional) Number of simulations to be run, 1 Million by default.
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

#' Calculates Type I Error Rate From Two-Sided Bounds
#'
#' Given bounds for a two sided test on uniform order statistics, this computes
#' the Type I Error Rate \eqn{\alpha} using a binary search.
#'
#' @param lower_bounds Numeric vector where the ith component is the lower bound
#' for the ith order statistic. The values must be in ascending order.
#' @param upper_bounds Numeric vector where the ith component is the lower bound
#' for the ith order statistic. The values must be in ascending order and
#' the ith component must be larger than the ith component of the lower bounds.
#'
#' @return Type I Error Rate \eqn{\alpha}
#'
#' @examples
#' get_level_from_bounds_two_sided(lower_bounds = c(.1, .5), upper_bounds = c(.6, .9))
#' 
#' @importFrom rlang .data
#'
#' @useDynLib qqconf jointlevel_twosided
#'
#' @export
get_level_from_bounds_two_sided <- function(lower_bounds,
                                           upper_bounds) {

  check_bounds_two_sided(lower_bounds, upper_bounds)
  n <- length(lower_bounds)
  b_vec <- c(lower_bounds, upper_bounds)
  h_g_df <- data.frame(b = b_vec, h_or_g = c(rep(0, n), rep(1, n)))
  h_g_df <- dplyr::arrange(h_g_df, .data$b)
  b_vec <- h_g_df$b
  bound_id <- h_g_df$h_or_g
  out <- 0.0
  # Below, a C routine is called for speed.
  res <- .C("jointlevel_twosided", b_vec = as.double(b_vec),
            bound_id = as.integer(bound_id), num_points = as.integer(n),
            out = as.double(out))
  return(res$out)

}

#' Calculates Local Bounds For Two-Sided Test
#'
#' Given the size of the dataset and the desired Type I Error rate \eqn{\alpha},
#' this calculates two-sided testing bounds for each order statistic.
#' These bounds are created using "equal local levels," which means that
#' each bound pair corresponds to a test that has a Type I Error rate of \eqn{\eta}
#' on its corresponding order statistic.
#'
#' @param alpha Desired Type I Error rate of the complete test.
#' @param n Size of the dataset.
#' @param tol (Optional) Relative tolerance of the \eqn{\alpha} level of the
#' simultaneous test. Defaults to 1e-8.
#' @param max_it (Optional) Maximum number of iterations of Binary Search Algorithm
#' used to find the bounds. Defaults to 100 which should be much larger than necessary
#' for a reasonable tolerance.
#' @param method (Optional) Parameter indicating if the calculation should be done using a highly
#' accurate approximation, "approximate", or if the calculations should be done using an exact
#' binary search calculation, "search". The approximate method is only usable for alpha values
#' of .1, .05, and .01.
#'
#' @return A list with components
#' \itemize{
#'   \item lower_bound - Numeric vector containing the lower bounds of the test of each order statistic.
#'   \item upper_bound - Numeric vector containing the upper bounds of the test of each order statistic.
#'   \item x - Expectation of each order statistic. This is the x-axis for the bounds if used in a qq-plot.
#'   \item local_level - Type I Error rate of each individual test on the order statistic. It is equal for all tests.
#' }
#'
#' @examples
#' get_bounds_two_sided(alpha = .05, n = 100, tol = 1e-6, max_it = 50)
#'
#'
#' @export
get_bounds_two_sided <- function(alpha,
                                n,
                                tol = 1e-8,
                                max_it = 100,
                                method=c("approximate", "search")) {

  `%>%` <- magrittr::`%>%`
  n_param <- n

  # Approximations are only available for alpha = .05 or alpha = .01
  if (method == "search") {

    eta_high <- -log(1 - alpha) / (2 * log(log(n)) * log(n)) # this is the asymptotic level of eta
    eta_low <- alpha / n
    eta_curr <- eta_low + (eta_high - eta_low) / 2
    n_it <- 0

    while (n_it < max_it) {

      n_it <- n_it + 1
      h_vals <- stats::qbeta(eta_curr / 2, 1:n, n:1)
      g_vals <- stats::qbeta(1 - (eta_curr / 2), 1:n, n:1)
      test_alpha <- 1 - get_level_from_bounds_two_sided(h_vals, g_vals)

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

    if (!(alpha %in% c(.1, .05, .01))) {

      stop("The approximate method is only configured for alpha = .1, .05, .01. Consider setting method='search'")

    }

    if (alpha == .01) {

      lookup_table <- alpha_01_df

    } else if (alpha == .05) {

      lookup_table <- alpha_05_df

    }

    if (alpha == .01 || alpha == .05) {

      if (n %in% lookup_table$n) {

        eta_df <- lookup_table %>%
          dplyr::filter(n == n_param)

        eta <- eta_df$local_level[1]

      } else if ((alpha == .05 && n > 5 * (10 ^ 4)) || (alpha == .01 && n > 10 ^ 5)) {

        eta <- get_asymptotic_approx_corrected_alpha(n, alpha)

      }
      else {

        # Do linear interpolation
        larger_n_df <- lookup_table %>%
          dplyr::filter(n > n_param) %>%
          dplyr::arrange(n) %>%
          dplyr::slice_head()

        larger_n <- larger_n_df$n[1]
        larger_n_eta <- larger_n_df$local_level[1]

        smaller_n_df <- lookup_table %>%
          dplyr::filter(n < n_param) %>%
          dplyr::arrange(n) %>%
          dplyr::slice_tail()

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

  alpha_vec <- seq(from = 1, to = n, by = 1)
  beta_vec <- n - alpha_vec + 1
  order_stats_mean <- alpha_vec / (alpha_vec + beta_vec)

  return(list(lower_bound = h_vals,
              upper_bound = g_vals,
              x = order_stats_mean,
              local_level = eta))

}

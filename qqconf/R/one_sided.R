#' Check Validity of One-Sided Bounds
#'
#' Given bounds for a one sided test, this checks that none of
#' the bounds fall outside of [0, 1].
#'
#' @param lower_bounds Numeric vector where the ith component is the lower bound
#' for the ith order statistic.
#' @param upper_bounds Numeric vector where the ith component is the lower bound
#' for the ith order statistic.
#'
#' @return None
#'
#' @examples
#' check_bounds_two_sided(lower_bounds = c(.1, .5, .8))
check_bounds_one_sided <- function(bounds) {
  
  if(any(bounds > 1) || any(bounds < 0)) {
    
    stop("Not all bounds between 0 and 1")
    
  }
  
}


#' Calculates Type I Error Rate From One-Sided Bounds
#'
#' Given bounds for a one sided test on uniform order statistics, this computes
#' the Type I Error Rate \eqn{\alpha} using a binary search. The function also
#' allows for an approximate computation for speed with some relative tolerance.
#'
#' @param bounds Numeric vector where the ith component is the lower bound
#' for the ith order statistic. The values must be in ascending order.
#'
#' @return Type I Error Rate \eqn{\alpha}
#'
#' @useDynLib qqconf jointlevel_onesided
#'
#' @examples
#' get_level_from_bounds_one_sided(bounds = c(.1, .5, .8))
#'
#' @export
get_level_from_bounds_one_sided <- function(bounds) {
  
  check_bounds_one_sided(bounds)
  n <- length(bounds)
  n_factorial_vec <- factorial(seq(0, n))
  bounds <- append(bounds, 1) # adding h_(n+1) = 1
  
  if (n == 1) {
    
    return(bounds[1])
    
  } else {
    
    c_0_1 <- (1 - bounds[2]) ^ n
    c_1_1 <- n * (bounds[2] - bounds[1]) * ((1 - bounds[2]) ^ (n - 1))
    c_prev <- c(c_0_1, c_1_1)
    
    for (k in seq(from = 2, to = (n - 1), by = 1)) {
      
      denom_vec <- (1 - bounds[k]) ^ seq((n - (k - 1)), n)
      num_vec <- (bounds[k + 1] - bounds[k]) ^ seq(0, (k - 1))
      c_curr <- rep(0, (k + 1)) 
      c_curr[1] <- (1 - bounds[k + 1]) ^ n
      
      for (j in seq(from = 1, to = (k), by = 1)) {
        
        jth_factorial_vec <- n_factorial_vec[(n + 1) : (n - j + 1)] / 
          (n_factorial_vec[(j + 1) : 1] * n_factorial_vec[(n - j + 1)])
        
        jth_quotient_vec <- ((1 - bounds[k + 1]) ^ (n - j)) * num_vec[(j + 1) : 1] / 
          denom_vec[length(denom_vec) : (length(denom_vec) - j)]
        
        c_curr[j + 1] <- sum(c_prev[1:(j)] * jth_factorial_vec * jth_quotient_vec)
        
      }
      
      c_prev <- c_curr
      
    }
    
    return(1 - sum(c_prev))
    
  }
  
}

#' Calculates Local Bounds For One-Sided Test
#'
#' Given the size of the dataset and the desired Type I Error rate \eqn{\alpha},
#' this calculates one-sided testing bounds for each order statistic.
#' These bounds are created using "equal local levels," which means that
#' each bound corresponds to a test that has a Type I Error rate of \eqn{\eta}
#' on its corresponding order statistic. Note that this code may be very slow for
#' values of n.
#'
#' @param alpha Desired Type I Error rate of the complete test.
#' @param n Size of the dataset.
#' @param tol (Optional) Relative tolerance of the \eqn{\alpha} level of the
#' simultaneous test. Defaults to 1e-8.
#' @param max_it (Optional) Maximum number of iterations of Binary Search Algorithm
#' used to find the bounds. Defaults to 100 which should be much larger than necessary
#' for a reasonable tolerance.
#'
#' @return A list with components
#' \itemize{
#'   \item bound - Numeric vector containing the bounds of the test of each order statistic.
#'   \item x - Expectation of each order statistic. This is the x-axis for the bounds if used in a qq-plot.
#'   \item local_level - Type I Error rate of each individual test on the order statistic. It is equal for all tests.
#' }
#'
#' @examples
#' get_level_from_bounds_onesided(alpha = .05, n = 100, tol = 1e-6, max_it = 50)
#'
#'
#' @export
get_bounds_onesided <- function(alpha, n, tol = 1e-8, max_it = 100) {
  
  eta_high <- -log(1 - alpha) / (2 * log(log(n)) * log(n))
  eta_low <- alpha / n
  eta_curr <- eta_low + (eta_high - eta_low) / 2
  n_it <- 0
  rel_err_tol <- (1 / alpha) * tol * (10 ^ (-4))
  bin_search_tol <- tol - rel_err_tol
  
  while (n_it < max_it) {
    
    n_it <- n_it + 1
    h_vals <- qbeta(eta_curr, 1:n, n:1)
    
    test_alpha <- 1 - get_level_from_bounds_one_sided(h_vals, approx = FALSE)
    
    if (abs(test_alpha - alpha) / alpha <= bin_search_tol) break
    
    if (test_alpha > alpha) {
      
      eta_high <- eta_curr
      eta_curr <- eta_curr - (eta_curr - eta_low) / 2
      
    } else if (test_alpha < alpha) {
      
      eta_low <- eta_curr
      eta_curr <- eta_curr + (eta_high - eta_curr) / 2
      
    }
    
  }
  
  alpha_vec <- seq(from = 1, to = n, by = 1)
  beta_vec <- n - alpha_vec + 1
  order_stats_mean <- alpha_vec / (alpha_vec + beta_vec)
  
  return(list(bound = h_vals, x = order_stats_mean, local_level = eta_curr))
  
}





jointlevel_onesided_exact_fast <- function(crit_vals) {
  
  n <- length(crit_vals)
  n_factorial_vec <- factorial(seq(0, n))
  crit_vals <- append(crit_vals, 1) # adding h_(n+1) = 1
  if (is.unsorted(crit_vals)) {
    
    stop("crit.vals must satisfy 0 < h_1 < h_2 < ... h_n < 1")
    
  }
  
  if (n == 1) {
    
    return(crit_vals[1])
    
  } else {
    
    c_0_1 <- (1 - crit_vals[2]) ^ n
    c_1_1 <- n * (crit_vals[2] - crit_vals[1]) * ((1 - crit_vals[2]) ^ (n - 1))
    c_prev <- c(c_0_1, c_1_1)
    
    for (k in seq(from = 2, to = (n - 1), by = 1)) {
      
      denom_vec <- (1 - crit_vals[k]) ^ seq((n - (k - 1)), n)
      num_vec <- (crit_vals[k + 1] - crit_vals[k]) ^ seq(0, (k - 1))
      c_curr <- rep(0, (k + 1)) 
      c_curr[1] <- (1 - crit_vals[k + 1]) ^ n
      
      # this was a k - 1
      for (j in seq(from = 1, to = (k), by = 1)) {
        
        jth_factorial_vec <- n_factorial_vec[(n + 1) : (n - j + 1)] / 
          (n_factorial_vec[(j + 1) : 1] * n_factorial_vec[(n - j + 1)])
        
        jth_quotient_vec <- ((1 - crit_vals[k + 1]) ^ (n - j)) * num_vec[(j + 1) : 1] / 
          denom_vec[length(denom_vec) : (length(denom_vec) - j)]
        
        c_curr[j + 1] <- sum(c_prev[1:(j)] * jth_factorial_vec * jth_quotient_vec)
        
      }
      
      c_prev <- c_curr
      
    }
    
    return(1 - sum(c_prev))
    
  }
  
}
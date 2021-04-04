#' Check Validity of One-Sided Bounds
#'
#' Given bounds for a one sided test, this checks that none of
#' the bounds fall outside of [0, 1].
#'
#' @param upper_bounds Numeric vector where the ith component is the upper bound
#' for the ith order statistic.
#'
#' @return None
#'
check_bounds_one_sided <- function(upper_bounds) {
  
  if(any(upper_bounds > 1) || any(upper_bounds < 0)) {
    
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
#' @examples
#' get_level_from_bounds_one_sided(bounds = c(.1, .5, .8))
#'
#' @export
get_level_from_bounds_one_sided <- function(bounds) {
  
  check_bounds_one_sided(bounds)
  n <- length(bounds)
  bounds <- append(bounds, 1) # adding h_(n+1) = 1
  if (is.unsorted(bounds)) {
    
    stop("crit.vals must satisfy 0 < h_1 < h_2 < ... h_n < 1")
    
  }
  
  if (n == 1) {
    
    return(bounds[1])
    
  } else {
    
    lgamma_vec <- numeric(n + 1)
    
    lgamma_vec[1] <- 0
    for (i in seq(2, n + 1)) {
      
      lgamma_vec[i] <- lgamma_vec[i - 1] + log(i - 1)
      
    }
    
    cc <- c((1 - bounds[2]) ^ n, exp(log(n) + log(bounds[2] - bounds[1]) + (n - 1) * log(1 - bounds[2])), 0)
    cn <- cc
    
    if (n > 2) for(k in c(2:(n-1))) {
      
      cn[1] <- (1 - bounds[k + 1]) ^ n
      
      for (j in c(1:(k-1))) {
        
        cn[j + 1] <- exp(log(cc[j + 1]) + (n - j) * log(1 - bounds[k + 1]) - (n - j) * log(1 - bounds[k]))
        
        for(l in c(0:(j-1))) {
          
          cn[j + 1] <- cn[j + 1] + exp(log(cc[l + 1]) + (j - l) * log(bounds[k + 1] - bounds[k]) + 
                          (n - j) * log(1 - bounds[k + 1]) - (n - l) * log(1 - bounds[k]) + lgamma_vec[n - l + 1] - lgamma_vec[j - l + 1] - lgamma_vec[n - j + 1])
          
        } 
      
      }
      
      cn[k + 1] <- 0
      for(l in c(0:(k-1))) {
        
        cn[k + 1] <- cn[k + 1] + exp(log(cc[l + 1]) + (k - l) * log(bounds[k + 1] - bounds[k]) + (n - k) * log(1 - bounds[k + 1]) - (n - l) *
                        log(1 - bounds[k]) + lgamma_vec[n - l + 1] - lgamma_vec[k - l + 1] - lgamma_vec[n - k + 1])
        
      } 
      
      cc <- cn
      cn <- c(cn, 0)
    }
    
    cn[n + 1] <- 0
    return(1-sum(cc))

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
#' get_bounds_one_sided(alpha = .05, n = 10, tol = 1e-6, max_it = 50)
#'
#'
#' @export
get_bounds_one_sided <- function(alpha, n, tol = 1e-6, max_it = 100) {
  
  if (alpha > .99 || n < 7) {
    
    eta_high <- alpha
    
  } else {
    
    eta_high <- -log(1 - alpha) / (2 * log(log(n)) * log(n))
    
  }
  
  eta_low <- alpha / n
  eta_curr <- eta_low + (eta_high - eta_low) / 2
  n_it <- 0
  
  while (n_it < max_it) {
    
    n_it <- n_it + 1
    h_vals <- stats::qbeta(eta_curr, 1:n, n:1)
    
    test_alpha <- get_level_from_bounds_one_sided(h_vals)
    
    if (abs(test_alpha - alpha) / alpha <= tol) break
    
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


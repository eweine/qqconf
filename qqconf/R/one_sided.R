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
  
  if(any(upper_bounds >= 1) || any(upper_bounds <= 0)) {
    
    stop("Not all bounds between 0 and 1 (exclusive)")
    
  }
  
  if(is.unsorted(upper_bounds)) {
    
    stop("Bounds are not sorted")
    
  }
  
  if(any(duplicated(upper_bounds))) {
    
    stop("Not all values of the bounds are distinct")
    
  }
  
}


#' Calculates Global Significance Level From Simultaneous One-Sided Bounds for Rejection Region
#'
#' For a one-sided test of uniformity of i.i.d. observations on the unit interval, 
#' this function will determine the significance level as a function of the rejection region. 
#' Suppose \eqn{n} observations are drawn i.i.d. from some CDF F(x) on the unit interval,
#' and it is desired to test the null hypothesis that F(x) = x for all x in (0, 1) against 
#' the one-sided alternative F(x) > x. Suppose the acceptance region for the test is 
#' described by a set of lower bounds, one for each order statistic.
#' Given the lower bounds, this function calculates the significance level of the test where the
#' null hypothesis is rejected if at least one of the order statistics 
#' falls below its corresponding lower bound.
#'
#' @param bounds Numeric vector where the ith component is the lower bound
#' for the ith order statistic. The components must be distinct values 
#' in (0, 1) that are in ascending order.
#'
#' @return Global significance level
#'
#' @examples
#' # For X1, X2, X3 i.i.d. unif(0, 1), 
#' # calculate 1 - P(X(1) > .1 and X(2) > .5 and X(3) > .8),
#' # where X(1), X(2), and X(3) are the order statistics.
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



#' Calculates Rejection Region of One-Sided Equal Local Levels Test
#'
#' The context is that n i.i.d. observations are assumed to be drawn
#' from some distribution on the unit interval with c.d.f. F(x), and it is 
#' desired to test the null hypothesis that F(x) = x for all x in (0,1), 
#' referred to as the "global null hypothesis," against the alternative F(x) > x for at least one x in (0, 1).
#' An "equal local levels" test is used, in which each of the n order statistics is
#' tested for significant deviation from its null distribution by a one-sided test 
#' with significance level \eqn{\eta}.  The global null hypothesis is rejected if at 
#' least one of the order statistic tests is rejected at level eta, where eta is 
#' chosen so that the significance level of the global test is alpha.  
#' Given the size of the dataset n and the desired global significance level alpha, 
#' this function calculates the local level eta and the acceptance/rejection regions for the test.
#' The result is a set of lower bounds, one for each order statistic.  
#' If at least one order statistic falls below the corresponding bound, 
#' the global test is rejected. Note that the code may be slow for n > 500.
#'
#'
#'
#' @param alpha Desired global significance level of the test.
#' @param n Size of the dataset.
#' @param tol (Optional) Relative tolerance of the \code{alpha} level of the
#' simultaneous test. Defaults to 1e-6.
#' @param max_it (Optional) Maximum number of iterations of Binary Search Algorithm
#' used to find the bounds. Defaults to 100 which should be much larger than necessary
#' for a reasonable tolerance.
#'
#' @return A list with components
#' \itemize{
#'   \item bound - Numeric vector of length \code{n} containing the lower bounds of the acceptance regions for the test of each order statistic.
#'   \item x - Numeric vector of length \code{n} containing the expectation of each order statistic. These are the x-coordinates for the bounds if used in a qq-plot.
#'   The value is \code{c(1:n) / (n + 1)}.
#'   \item local_level - Significance level \eqn{\eta} of the local test on each individual order statistic. It is equal for all order
#'   statistics and will be less than \code{alpha} for all \code{n} > 1.
#' }
#'
#' @examples
#' get_bounds_one_sided(alpha = .05, n = 10, max_it = 50)
#'
#'
#' @export
get_bounds_one_sided <- function(alpha, n, tol = 1e-6, max_it = 100) {
  
  if (alpha >= 1 || alpha <= 0) {
    
    stop("alpha must be between 0 and 1 (exclusive).")
    
  }
  
  if (as.integer(n) != n) {
    
    stop("n must be an integer")
    
  }
  
  if (n < 1) {
    
    stop("n must be greater than 0")
    
  }
  
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


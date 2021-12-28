#' Estimate parameters of distribution given set of observations
#'
#'   For the normal distribution, we estimate the mean as the median and the standard deviation as \eqn{Sn} from the paper by Rousseeuw and Croux 1993
#'   "Alternatives to the Median Absolute Deviation". For all other distributions besides normal,
#'   the code uses MLE to estimate the parameters.
#'
#' @param obs The observed data.
#' @param distribution The quantile  or probability function for the specified distribution given as a string.
#' For example \code{"qnorm"} or \code{"pnorm"}.
#'
#' @return list of estimated parameters
estimate_dparams <- function(distribution, obs) {

  first_letter <- substr(distribution, 1, 1)

  if (!(first_letter %in% c("q", "p"))) {

    stop("distribution is not recognized. Please pass in a distribution function that begins with a 'q' or 'p'.
          Assigning a distribution function to another variable will not work, please pass in the distribution
          function itself")

  }

  dist_name <- substr(distribution, 2, nchar(distribution))

  # equivalence between base R and MASS::fitdistr distribution names
  corresp <- function(distributionName) {
    switch(
      distributionName,
      beta = "beta",
      cauchy = "cauchy",
      chisq = "chi-squared",
      exp = "exponential",
      f = "f",
      gamma = "gamma",
      geom = "geometric",
      lnorm = "log-normal",
      logis = "logistic",
      norm = "normal",
      nbinom = "negative binomial",
      pois = "poisson",
      t = "t",
      weibull = "weibull",
      NULL
    )
  }

  # initial value for some distributions
  initVal <- function(distributionName) {
    switch(
      distributionName,
      beta = list(shape1 = 1, shape2 = 1),
      chisq = list(df = 1),
      f = list(df1 = 1, df2 = 2),
      t = list(df = 1),
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

  return(dparams)

}

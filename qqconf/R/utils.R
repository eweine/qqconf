#' Shorthand for two numerical comparisons
#'
#' @param x numeric value
#' @param gte lower bound
#' @param lte upper bound
#'
#' @return boolean
between <- function(x, gte, lte) {

  x >= gte & x <= lte

}

#' Convert R Distribution Function to MASS Distribution Name
#'
#' @param distr R distribution function (e.g. qnorm or pnorm)
#' @param band_type one of "qq" (for quantile functions) or "pp" (
#' for probability functions).
#'
#' @return string of MASS distribution name
get_mass_name_from_distr <- function(distr, band_type) {

  if (band_type == "qq") {

    corresp_dn <- function(distributionName) {
      switch(
        distributionName,
        qbeta = "beta",
        qcauchy = "cauchy",
        qchisq = "chi-squared",
        qexp = "exponential",
        qf = "f",
        qgamma = "gamma",
        qgeom = "geometric",
        qnorm = "normal",
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

  } else if (band_type == "pp") {

    corresp_dn <- function(distributionName) {
      switch(
        distributionName,
        pbeta = "beta",
        pcauchy = "cauchy",
        pchisq = "chi-squared",
        pexp = "exponential",
        pf = "f",
        pgamma = "gamma",
        pgeom = "geometric",
        pnorm = "normal",
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

  }

  dist <- corresp_dn(distr)
  if(is.null(dist)) {

    stop("Unknown distribution provided")

  }
  return(dist)

}

#' Estimate Parameters from Data
#'
#' For select distributions, parameters are estimated from data. Generally,
#' the MLEs are used. However, for the normal distribution we use robust
#' estimators.
#'
#' @param distr_name \code{MASS} name of distribution
#' @param obs observation vector
#'
#' @importFrom stats sd
#'
#' @return list of distribution parameters
estimate_params_from_data <- function(distr_name, obs) {

  # initial value for some distributions
  initVal <- function(distributionName) {
    switch(
      distributionName,
      "beta" = list(shape1 = 1, shape2 = 1),
      "chi-squared" = list(df = 1),
      "f" = list(df1 = 1, df2 = 2),
      NULL
    )
  }

    if(is.null(initVal(distr_name))) {
      if(distr_name == "normal") {

        # use special estimators for the normal distribution
        dparams <- list()
        dparams['mean'] <- median(x = obs)

        if (requireNamespace("robustbase", quietly = TRUE)) {

          dparams['sd'] <- robustbase::Sn(x = obs)

        } else {

          warning(
            "Using an MLE to estimate the sd of a normal. Resulting bands may be highly conservative.
            To use a more robust estimator, please install robustbase and rerun the code",
            call. = FALSE,
            immediate. = TRUE
            )
          dparams['sd'] <- sd(obs)

        }

      } else {

        suppressWarnings({
          dparams <- MASS::fitdistr(x = obs, densfun = distr_name)$estimate
        })

      }

    } else {

      suppressWarnings({
        dparams <- MASS::fitdistr(
          x = obs,
          densfun = distr_name,
          start = initVal(distr_name))$estimate
      })

    }


  return(dparams)

}

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
        dparams['sd'] <- robustbase::Sn(x = obs)

      } else if (distr_name == "t") {

        dparams <- list()
        suppressWarnings({
          param_ests <- MASS::fitdistr(
            x = obs,
            densfun = "t",
            start = list(df = 2, m = 1),
            s = 1
          )$estimate
        })

        dparams['ncp'] <- param_ests['m']
        dparams['df'] <- param_ests['df']

      }
      else {

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

#' Get Best Available Method for Probability Points
#'
#' Determines name of best method for obtaining expected points for a qq or pp
#' plot.
#'
#' @param dist_name character name of distribution
#'
#' @return character name of best expected points method
#'
get_best_available_prob_pts_method <- function(dist_name) {

  if (dist_name %in% c("qunif", "punif")) {

    return("uniform")

  } else if (dist_name %in% c("qnorm", "pnorm")) {

    return("normal")

  } else {

    return("median")

  }

}

#' Get Quantile for First and Last Point of QQ or PP Plot
#'
#' @param exp_pts_method method used to derive expected points
#' @param n sample size
#'
#' @return list with low and high point
#'
get_extended_quantile <- function(exp_pts_method, n) {

  if (exp_pts_method == "uniform") {

    high_pt <- 1 - (1 / max(n * 1.25, n + 2))
    low_pt <- 1 / max(n * 1.25, n + 2)

  } else if (exp_pts_method == "median"){

    new_samp_size <- floor(n * 1.02)
    ppoints_adj <- ppoints(new_samp_size)
    low_pt <- ppoints_adj[1]
    high_pt <- ppoints_adj[new_samp_size]

  } else if (exp_pts_method == "normal") {

    new_samp_size <- floor(n * 1.3)
    ppoints_adj <- ppoints(new_samp_size)
    low_pt <- ppoints_adj[1]
    high_pt <- ppoints_adj[new_samp_size]

  }

  return(
    list(
      low_pt = low_pt,
      high_pt = high_pt
    )
  )

}

#' Get Quantile Distribution from Probability Distribution
#'
#' @param dname probability distribution (e.g. \code{pnorm})
#'
#' @return quantile distribution (e.g. \code{qnorm}).
#'
get_qq_distribution_from_pp_distribution <- function(dname) {

  dname_no_p <- substr(dname, 2, nchar(dname))
  eval(parse(text = paste0("q", dname_no_p)))

}

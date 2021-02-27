#' Probability-probability testing bands
#'
#' Draws probability-probability testing bands, with an additional difference option
#'
#' @import ggplot2
#'
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_ribbon
#'
#' @param distribution The quantile function for the specified distribution. Defaults to pnorm.
#' Custom distributions are allowed so long as all parameters are supplied in dparams.
#' @param dparams List of additional parameters for the quantile function of the distribution
#'   (e.g. df=1). Will be estimated if not provided and an appropriate estimation procedure exists.
#'   For the normal distribution, we estimate the mean as the median and the standard deviation as \eqn{Sn} from the paper by Rousseeuw and Croux 1993
#'   "Alternatives to the Median Absolute Deviation". For all other distributions,
#'   the code uses MLE to estimate the parameters. Note that estimation is not implemented for custom distributions, so all
#'   parameters of the distribution must be provided by the user.
#' @param difference Whether to plot the difference between the observed and
#'   expected values on the vertical axis.
#' @param method Method for simultaneous testing bands. Must be either "ell", which applies a level \eqn{\eta} pointwise
#' test to each order statistic such that the Type I error of the global test is \eqn{\alpha}, or "ks" to apply a 
#' Kolmogorov-Smirnov test. For \eqn{\alpha} = .01, .05, and .1, "ell" is recommended.
#' @param bounds_params List of optional parameters for get_bounds_two_sided
#'   (i.e. tol, max_it, method).
#' @param alpha Type I error of global test of if the data comes from the reference distribution.
#'
#' @examples
#' # generate random Normal data
#' set.seed(0)
#' smp <- data.frame(norm = rnorm(100), exp = rexp(100))
#'
#' # Normal P-P plot of Normal data
#' gg <- ggplot(data = smp, mapping = aes(sample = norm)) +
#'  stat_pp_band() +
#'  stat_pp_line() +
#'  stat_pp_point() +
#'  labs(x = "Probability Points", y = "Cumulative Probability")
#' gg
#'
#' # Shifted Normal P-P plot of Normal data
#' dp <- list(mean = 1.5)
#' gg <- ggplot(data = smp, mapping = aes(sample = norm)) +
#'  stat_pp_band(dparams = dp) +
#'  stat_pp_line() +
#'  stat_pp_point(dparams = dp) +
#'  labs(x = "Probability Points", y = "Cumulative Probability")
#' gg
#'
#' # Exponential P-P plot of Exponential data
#' di <- "exp"
#' gg <- ggplot(data = smp, mapping = aes(sample = exp)) +
#'  stat_pp_band(distribution = di) +
#'  stat_pp_line() +
#'  stat_pp_point(distribution = di) +
#'  labs(x = "Probability Points", y = "Cumulative Probability")
#' gg
#'
#' # Normal P-P plot of mean ozone levels (airquality dataset)
#' dp <- list(mean = 38, sd = 27)
#' gg <- ggplot(data = airquality, mapping = aes(sample = Ozone)) +
#'  stat_pp_band(dparams = dp) +
#'  stat_pp_line() +
#' 	stat_pp_point(dparams = dp) +
#'  labs(x = "Probability Points", y = "Cumulative Probability")
#' gg
#'
#' @export
stat_pp_band <- function(
	mapping = NULL,
	data = NULL,
	geom = "ribbon",
	position = "identity",
	na.rm = TRUE,
	show.legend = NA,
	inherit.aes = TRUE,
	distribution = "pnorm",
	dparams = list(),
	method = "ell",
	alpha = .05,
	difference = FALSE,
	...
) {
  
  conf <- 1 - alpha
	# error handling
	if (!(distribution %in% c(
		"pbeta",
		"pcauchy",
		"pchisq",
		"pexp",
		"pf",
		"pgamma",
		"pgeom",
		"plnorm",
		"plogis",
		"pnorm",
		"pnbinom",
		"ppois",
		"pt",
		"pweibull")) &
		length(dparams) == 0 &
		table(sapply(formals(eval(parse(text = distribution))), typeof))["symbol"] > 1) {
		stop(
			"MLE is currently not supported for custom distributions.\n",
			"Please provide all the custom distribution parameters to 'dparams'.",
			call. = FALSE
		)
	}
	if (conf < 0 | conf > 1) {
		stop("Please provide a valid alpha level for the bands: ",
				 "'conf' must be between 0 and 1.",
				 call. = FALSE)
	}

	# vector with common discrete distributions
	discreteDist <- c("binom", "geom", "hyper", "multinom", "nbinom", "pois")

	if (distribution %in% discreteDist) geom <- "errorbar"

	ggplot2::layer(
		mapping = mapping,
		stat = StatPpBand,
		geom = geom,
		position = position,
		show.legend = show.legend,
		inherit.aes = inherit.aes,
		params = list(
			na.rm = na.rm,
			distribution = distribution,
			dparams = dparams,
		  method = method,
			alpha = alpha,
			discrete = distribution %in% discreteDist,
			difference = difference,
			...
		)
	)
}

#' StatPpBand
#'
#' @keywords internal
#' @usage NULL
#' @export
StatPpBand <- ggplot2::ggproto(
	`_class` = "StatPpBand",
	`_inherit` = ggplot2::Stat,

	default_aes = ggplot2::aes(
		x = ..x..,
		ymin = ..lower..,
		ymax = ..upper..
	),

	required_aes = c("sample"),

	compute_group = {
		function(data,
						 self,
						 scales,
						 distribution,
						 dparams,
						 bandType,
						 conf,
						 discrete,
						 difference) {
			# distributional functions
			pFunc <- distribution

			smp <- data$sample
			n <- length(smp)
			probs <- ppoints(n)

			# automatically estimate parameters with MLE, only if no parameters are
			# provided with dparams and there are at least one distributional parameter
			# without a default value
			if(length(dparams) == 0) {
				# equivalence between base R and MASS::fitdistr distribution names
				corresp <- function(distName) {
					switch(
						distName,
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
						t = dt,
						weibull = "weibull",
						NULL
					)
				}

				# initial value for some distributions
				initVal <- function(distName) {
					switch(
						distName,
						beta = list(shape1 = 1, shape2 = 1),
						chisq = list(df = 1),
						f = list(df1 = 1, df2 = 2),
						t = list(df = 1),
						NULL
					)
				}

				suppressWarnings({
					if(!is.null(corresp(distribution))) {
						if(is.null(initVal(distribution))) {
							dparams <- MASS::fitdistr(x = smp, densfun = corresp(distribution))$estimate
						} else {
							dparams <- MASS::fitdistr(x = smp, densfun = corresp(distribution), start = initVal(distribution))$estimate
						}
					}
				})
			}

			# pointwise confidence bands based on normal confidence intervals
			if (method == "pointwise") {
			  probs <- ppoints(n)
			  stdErr <- (slope / do.call(pFunc, c(list(x = theoretical), dparams))) * sqrt(probs * (1 - probs) / n)
			  zCrit <- qnorm(p = (1 - (1 - conf) / 2))
			  
			  upper <- fittedValues + (stdErr * zCrit)
			  lower <- fittedValues - (stdErr * zCrit)
			}
			
			# using the DKW inequality for simultaneous bands
			if (method == "ks") {
			  probs <- ppoints(n)
			  epsilon <- sqrt((1 / (2 * n)) * log(2/(1-conf)))
			  lp <- pmax(probs - epsilon, rep(0, n))
			  up <- pmin(probs + epsilon, rep(1, n))
			  lower <- intercept + slope * do.call(pFunc, c(list(p = lp), dparams))
			  upper <- intercept + slope * do.call(pFunc, c(list(p = up), dparams))
			}
			
			# tail-sensitive confidence bands
			if (method == "ell") {
			  
			  global.bounds <- do.call(get_bounds_two_sided,
			                           c(list(alpha = 1 - conf, n = n), bounds_params))
			  lower <- do.call(pFunc, c(list(p = global.bounds$lower_bound), dparams))
			  upper <- do.call(pFunc, c(list(p = global.bounds$upper_bound), dparams))
			  
			}
			

			out <- data.frame(
				x = probs,
				upper = upper,
				lower = lower,
				fill = if (is.null(data$fill)) rgb(.6, .6, .6, .5) else data$fill
			)

			if (discrete) {
				out$colour <- rgb(.5, .5, .5)
				# create a data.frame with unique rows
				out <- dplyr::summarize(
					dplyr::group_by(out, x, fill, colour),
					upper = max(upper),
					lower = min(lower)
				)
				out <- as.data.frame(out)
			}

			# difference the confidence bands by keeping the same distance from the
			# identity line, which now is a line centered on y = 0
			if (difference) {
				out$upper <- out$upper - probs
				out$lower <- out$lower - probs
			}

			out
		}
	}
)

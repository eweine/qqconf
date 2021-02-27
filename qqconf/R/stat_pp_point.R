#' Probability-probability points
#'
#' Draws probability-probability points.
#'
#' @import ggplot2
#'
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_point
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
#'
#'
#' @examples
#' # generate random Normal data
#' set.seed(0)
#' smp <- data.frame(norm = rnorm(100))
#'
#' # Normal P-P plot of Normal data
#' gg <- ggplot(data = smp, mapping = aes(sample = norm)) +
#'  stat_pp_point() +
#'  labs(x = "Probability Points", y = "Cumulative Probability")
#' gg
#'
#' # Shifted Normal P-P plot of Normal data
#' dp <- list(mean = 1.5)
#' gg <- ggplot(data = smp, mapping = aes(sample = norm)) +
#'  stat_pp_point(dparams = dp) +
#'  labs(x = "Probability Points", y = "Cumulative Probability")
#' gg
#'
#' # Normal P-P plot of mean ozone levels (airquality dataset)
#' dp <- list(mean = 38, sd = 27)
#' gg <- ggplot(data = airquality, mapping = aes(sample = Ozone)) +
#' 	stat_pp_point(dparams = dp) +
#'  labs(x = "Probability Points", y = "Cumulative Probability")
#' gg
#'
#' @export
stat_pp_point <- function(
	mapping = NULL,
	data = NULL,
	geom = "point",
	position = "identity",
	na.rm = TRUE,
	show.legend = NA,
	inherit.aes = TRUE,
	distribution = "pnorm",
	dparams = list(),
	difference = FALSE,
	...
) {
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

	ggplot2::layer(
		mapping = mapping,
		stat = StatPpPoint,
		geom = geom,
		position = position,
		show.legend = show.legend,
		inherit.aes = inherit.aes,
		params = list(
			na.rm = na.rm,
			distribution = distribution,
			dparams = dparams,
			difference = difference,
			...
		)
	)
}

#' StatPpPoint
#'
#' @keywords internal
#' @usage NULL
#' @export
StatPpPoint <- ggplot2::ggproto(
	`_class` = "StatPpPoint",
	`_inherit` = ggplot2::Stat,

	default_aes = ggplot2::aes(x = ..theoretical.., y = ..sample..),

	required_aes = c("sample"),

	compute_group = function(data,
													 self,
													 scales,
													 distribution,
													 dparams,
													 difference) {
		# cumulative distributional function
		pFunc <- distribution

		smp <- sort(data$sample)
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

		# evaluate the cdf on the observed quantiles
		y <- do.call(pFunc, c(list(q = smp), dparams))

		if (difference) {
			# calculate new ys for the differenceed sample using the identity line
			dY <- y - probs

			out <- data.frame(sample = dY, theoretical = probs)
		} else {
			out <- data.frame(sample = y, theoretical = probs)
		}

		out
	}
)

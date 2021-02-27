#' Quantile-quantile lines
#'
#' Draws a quantile-quantile line, with an additional difference option.
#'
#'
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_path
#'
#' @param distribution The quantile function for the specified distribution. Defaults to qnorm.
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
#' # Normal Q-Q plot of Normal data
#' gg <- ggplot(data = smp, mapping = aes(sample = norm)) +
#'  stat_qq_line() +
#'  stat_qq_point() +
#'  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
#' gg
#'
#' # Exponential Q-Q plot of mean ozone levels (airquality dataset)
#' di <- "exp"
#' dp <- list(rate = 1)
#' gg <- ggplot(data = airquality, mapping = aes(sample = Ozone)) +
#'  stat_qq_line(distribution = di, dparams = dp) +
#'  stat_qq_point(distribution = di, dparams = dp) +
#'  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
#' gg
#'
#' # differenceed Exponential Q-Q plot of mean ozone levels
#' di <- "exp"
#' dp <- list(rate = 1)
#' de <- TRUE
#' gg <- ggplot(data = airquality, mapping = aes(sample = Ozone)) +
#'  stat_qq_line(distribution = di, difference = de) +
#'  stat_qq_point(distribution = di, difference = de) +
#'  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
#' gg
#'
#' @export
stat_qq_line <- function(
	mapping = NULL,
	data = NULL,
	geom = "path",
	na.rm = TRUE,
	position = "identity",
	show.legend = NA,
	inherit.aes = TRUE,
	distribution = "norm",
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
		table(sapply(formals(eval(parse(text = paste0("q", distribution)))), typeof))["symbol"] > 1) {
		stop(
			"MLE is currently not supported for custom distributions.\n",
			"Please provide all the custom distribution parameters to 'dparams'.",
			call. = FALSE
		)
	}

	ggplot2::layer(
		data = data,
		mapping = mapping,
		stat = StatQqLine,
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

#' StatQqLine
#'
#' @keywords internal
#' @usage NULL
#' @export
StatQqLine <- ggplot2::ggproto(
	`_class` = "StatQqLine",
	`_inherit` = ggplot2::Stat,

	required_aes = c("sample"),

	default_aes = ggplot2::aes(
		x = ..xline..,
		y = ..yline..
	),

	compute_group = {
		function(data,
						 self,
						 scales,
						 distribution,
						 dparams,
						 difference) {
			# distributional function
			qFunc <- eval(parse(text = paste0("q", distribution)))

			smp <- sort(data$sample)
			n <- length(smp)
			quantiles <- ppoints(n)

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

			theoretical <- do.call(qFunc, c(list(p = quantiles), dparams))

			if (difference) {
				out <- data.frame(xline = c(min(theoretical), max(theoretical)))
				out$yline <- 0
			} else {

				slope <- 1
				intercept <- 0

				out <- data.frame(xline = c(min(theoretical), max(theoretical)))
				out$yline <- slope * out$xline + intercept
			}

			out$size <- .8
			# out$colour <- if (is.null(data$colour)) rgb(.3, .3, .3) else rep(data$colour[1], 2)

			out
		}
	}
)

#' Quantile-quantile points
#'
#' Draws quantile-quantile points, with an additional difference option.
#'
#'
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_point
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
#' # Normal Q-Q plot of simulated Normal data
#' gg <- ggplot(data = smp, mapping = aes(sample = norm)) +
#'  stat_qq_point() +
#'  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
#' gg
#'
#' # Exponential Q-Q plot of mean ozone levels (airquality dataset)
#' di <- "exp"
#' dp <- list(rate = 1)
#' gg <- ggplot(data = airquality, mapping = aes(sample = Ozone)) +
#'  stat_qq_point(distribution = di, dparams = dp) +
#'  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
#' gg
#'
#' @export
stat_qq_point <- function(
	mapping = NULL,
	data = NULL,
	geom = "point",
	position = "identity",
	na.rm = TRUE,
	show.legend = NA,
	inherit.aes = TRUE,
	distribution = "norm",
	dparams = list(),
	difference = FALSE,
	...
) {
	# error handling
	if (!(distribution %in% c(
		"qbeta",
		"qcauchy",
		"qchisq",
		"qexp",
		"qf",
		"qgamma",
		"qgeom",
		"qlnorm",
		"qlogis",
		"qnorm",
		"qnbinom",
		"qpois",
		"qt",
		"qweibull")) &
		length(dparams) == 0 &
		table(sapply(formals(eval(parse(text = distribution))), typeof))["symbol"] > 1) {
		stop(
			"MLE is currently not supported for custom distributions.\n",
			"Please provide all the custom distribution parameters to 'dparams'.",
			call. = FALSE
		)
	}


	ggplot2::layer(
		data = data,
		mapping = mapping,
		stat = StatQqPoint,
		position = position,
		geom = geom,
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

#' StatQqPoint
#'
#' @keywords internal
#' @usage NULL
#' @export
StatQqPoint <- ggplot2::ggproto(
	`_class` = "StatQqPoint",
	`_inherit` = ggplot2::Stat,

	default_aes = ggplot2::aes(
		x = ..theoretical..,
		y = ..sample..
	),

	required_aes = c("sample"),

	compute_group = function(data,
													 self,
													 scales,
													 distribution,
													 dparams,
													 difference) {

		oidx <- order(data$sample)
		smp <- data$sample[oidx]
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
					qbeta = "beta",
					qcauchy = "cauchy",
					qchisq = "chi-squared",
					qexp = "exponential",
					qf = "f",
					qgamma = "gamma",
					qgeom = "geometric",
					qlnorm = "log-normal",
					qlogis = "logistic",
					qnorm = "normal",
					qnbinom = "negative binomial",
					qpois = "poisson",
					qt = dt,
					qweibull = "weibull",
					NULL
				)
			}

			# initial value for some distributions
			initVal <- function(distName) {
				switch(
					distName,
					qbeta = list(shape1 = 1, shape2 = 1),
					qchisq = list(df = 1),
					qf = list(df1 = 1, df2 = 2),
					qt = list(df = 1),
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

		theoretical <- do.call(distribution, c(list(p = quantiles), dparams))

		if (difference) {
			slope <- 1
			intercept <- 0

			# calculate new ys for the differenceed sample
			dSmp <- NULL
			for (i in 1:n) {
				lSmp <- slope * theoretical[i] + intercept
				dSmp[i] <- smp[i] - lSmp
			}

			smp <- dSmp
		}
		
		out <- data.frame(sample = smp, theoretical = theoretical)

		out
	}
)


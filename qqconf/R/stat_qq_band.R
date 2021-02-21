#' Quantile-quantile testing bands
#'
#' Draws quantile-quantile confidence bands, with an additional difference option.
#'
#' If any of the points of the qq-plot fall outside the simultaneous acceptance region for the selected
#' level alpha test, that means that we can reject the null hypothesis that the data are i.i.d. draws from the
#' specified distribution. If 'difference' is set to TRUE, the vertical axis plots the 
#' observed quantile minus expected quantile. Set pw.lty to a non-zero line type to plot
#' the pointwise bounds. If pointwise bands are used, then on average, alpha * n of the points will fall outside
#' the bounds under the null hypothesis, so the chance that the qq-plot has any points falling outside of the pointwise bounds
#' is typically much higher than alpha under the null hypothesis. For this reason, a simultaneous region is preferred. 
#'
#' @include stat_qq_line.R
#'
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_ribbon
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
#' @param method Method for simultaneous testing bands. Must be either "ell", which applies a level \eqn{\eta} pointwise
#' test to each order statistic such that the Type I error of the global test is \eqn{\alpha}, or "ks" to apply a 
#' Kolmogorov-Smirnov test. For \eqn{\alpha} = .01, .05, and .1, "ell" is recommended.
#' @param bounds_params List of optional parameters for get_bounds_two_sided
#'   (i.e. tol, max_it, method).
#' @param alpha Type I error of global test of if the data comes from the reference distribution.
#'
#' @export
stat_qq_band <- function(
	mapping = NULL,
	data = NULL,
	geom = "qq_band",
	position = "identity",
	na.rm = TRUE,
	show.legend = NA,
	inherit.aes = TRUE,
	distribution = "norm",
	dparams = list(),
	difference = FALSE,
	method = "ell",
	bounds_params = NULL,
	alpha = .05,
	...
) {
  
  conf <- 1 - alpha
	# error handling
	if (!(distribution %in% c(
		"beta",
		"cauchy",
		"chisq",
		"exp",
		"f",
		"gamma",
		"geom",
		"lnorm",
		"logis",
		"norm",
		"nbinom",
		"pois",
		"t",
		"weibull")) &
		length(dparams) == 0 &
		table(sapply(formals(eval(parse(text = paste0("q", distribution)))), typeof))["symbol"] > 1) {
		stop(
			"MLE is currently not supported for custom distributions.\n",
			"Please provide all the custom distribution parameters to 'dparams'.",
			call. = FALSE
		)
	}
	if (conf < 0 | conf > 1) {
		stop("Please provide a valid alpha value for the bands: ",
				 "'alpha' must be between 0 and 1.",
				 call. = FALSE)
	}
	method <- match.arg(method, c("pointwise", "ks", "ell"))

	# vector with common discrete distributions
	discreteDist <- c("binom", "geom", "nbinom", "pois")

	if (distribution %in% discreteDist) geom <- "errorbar"

	ggplot2::layer(
		data = data,
		mapping = mapping,
		stat = StatQqBand,
		geom = geom,
		position = position,
		show.legend = show.legend,
		inherit.aes = inherit.aes,
		params = list(
			na.rm = na.rm,
			distribution = distribution,
			dparams = dparams,
			difference = difference,
			method = method,
			bounds_params = bounds_params,
			conf = conf,
			discrete = distribution %in% discreteDist,
			...
		)
	)
}

#' StatQqBand
#'
#' @keywords internal
#' @usage NULL
#' @export
StatQqBand <- ggplot2::ggproto(
	`_class` = "StatQqBand",
	`_inherit` = StatQqLine,

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
						 difference,
						 method,
						 bounds_params,
						 conf,
						 discrete) {
			# distributional functions
			dFunc <- eval(parse(text = paste0("d", distribution)))
			qFunc <- eval(parse(text = paste0("q", distribution)))
			rFunc <- eval(parse(text = paste0("r", distribution)))

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
						t = "t",
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
						  if(corresp(distribution) == "normal") {
						    
						    # Use special estimators for the normal distribution
						    dparams <- c()
						    dparams['mean'] <- median(x = smp)
						    dparams['sd'] <- robustbase::Sn(x = smp)
						    
						  }
							dparams <- MASS::fitdistr(x = smp, densfun = corresp(distribution))$estimate
						} else {
							dparams <- MASS::fitdistr(x = smp, densfun = corresp(distribution), start = initVal(distribution))$estimate
						}
					}
				})
			}

			theoretical <- do.call(qFunc, c(list(p = quantiles), dparams))

			# inherit from StatQqLine
			xline <- self$super()$compute_group(data = data,
																					distribution = distribution,
																					dparams = dparams,
																					difference = FALSE)$xline
			yline <- self$super()$compute_group(data = data,
																					distribution = distribution,
																					dparams = dparams,
																					difference = FALSE)$yline

			slope <- diff(yline) / diff(xline)
			intercept <- yline[1L] - slope * xline[1L]

			fittedValues <- (slope * theoretical) + intercept

			# pointwise confidence bands based on normal confidence intervals
			if (method == "pointwise") {
				probs <- ppoints(n)
				stdErr <- (slope / do.call(dFunc, c(list(x = theoretical), dparams))) * sqrt(probs * (1 - probs) / n)
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
				lower <- intercept + slope * do.call(qFunc, c(list(p = lp), dparams))
				upper <- intercept + slope * do.call(qFunc, c(list(p = up), dparams))
			}

			# tail-sensitive confidence bands
			if (method == "ell") {
				
				global.bounds <- do.call(get_bounds_two_sided,
				                         c(list(alpha = 1 - conf, n = n), bounds_params))
				lower <- do.call(qFunc, c(list(p = global.bounds$lower_bound), dparams))
				upper <- do.call(qFunc, c(list(p = global.bounds$upper_bound), dparams))
				
			}

			out <- data.frame(
				x = theoretical,
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
			# stat_qq_line, which now is a line centered on y = 0
			if (difference) {
				out$upper <- out$upper - fittedValues
				out$lower <- out$lower - fittedValues
			}

			out
		}
	}
)

#' Quantile-quantile confidence bands
#'
#' Draws quantile-quantile confidence bands, with an additional difference option.
#'
#'
#' @include stat_qq_line.R
#'
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_ribbon
#'
#' @param distribution Character. Theoretical probability distribution function
#'   to use. Do not provide the full distribution function name (e.g.,
#'   \code{"dnorm"}). Instead, just provide its shortened name (e.g.,
#'   \code{"norm"}). If you wish to provide a custom distribution, you may do so
#'   by first creating the density, quantile, and random functions following the
#'   standard nomenclature from the \code{stats} package (i.e., for
#'   \code{"custom"}, create the \code{dcustom}, \code{pcustom},
#'   \code{qcustom}, and \code{rcustom} functions).
#' @param dparams List of additional parameters passed on to the previously
#'   chosen \code{distribution} function. If an empty list is provided (default)
#'   then the distributional parameters are estimated via MLE. MLE for custom
#'   distributions is currently not supported, so you must provide the
#'   appropriate \code{dparams} in that case.
#' @param difference Logical. Should the plot objects be differenceed? If \code{TRUE},
#'   the objects will be differenceed according to the reference Q-Q line. This
#'   procedure was described by Thode (2002), and may help reducing visual bias
#'   caused by the orthogonal distances from Q-Q points to the reference line.
#' @param bandType Character. Either \code{"pointwise"}, \code{"boot"}, \code{"ks"} or
#'   \code{"ts"}. \code{"pointwise"} constructs pointwise confidence bands based
#'   on Normal confidence intervals. \code{"boot"} creates pointwise confidence
#'   bands based on a parametric bootstrap; parameters are estimated with MLEs.
#'   \code{"ks"} constructs simultaneous confidence bands based on the Kolmogorov-Smirnov
#'   test. Finally, \code{"ts"} constructs tail-sensitive confidence bands, as
#'   described by Aldor-Noiman et al. (2013) (also, see 'Note' for
#'   limitations).
#' @param bounds.params List of optional parameters for get_bounds_two_sided
#'   (i.e. tol, max_it, method).
#' @param conf Numerical. Confidence level of the bands.
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
	bandType = "pointwise",
	bounds.params = NULL,
	conf = .95,
	...
) {
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
		stop("Please provide a valid confidence level for the bands: ",
				 "'conf' must be between 0 and 1.",
				 call. = FALSE)
	}
	bandType <- match.arg(bandType, c("pointwise", "ks", "equal_local_levels"))

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
			bandType = bandType,
			bounds.params = bounds.params,
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
						 bandType,
						 bounds.params,
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
			if (bandType == "pointwise") {
				probs <- ppoints(n)
				stdErr <- (slope / do.call(dFunc, c(list(x = theoretical), dparams))) * sqrt(probs * (1 - probs) / n)
				zCrit <- qnorm(p = (1 - (1 - conf) / 2))

				upper <- fittedValues + (stdErr * zCrit)
				lower <- fittedValues - (stdErr * zCrit)
			}

			# using the DKW inequality for simultaneous bands
			if (bandType == "ks") {
				probs <- ppoints(n)
				epsilon <- sqrt((1 / (2 * n)) * log(2/(1-conf)))
				lp <- pmax(probs - epsilon, rep(0, n))
				up <- pmin(probs + epsilon, rep(1, n))
				lower <- intercept + slope * do.call(qFunc, c(list(p = lp), dparams))
				upper <- intercept + slope * do.call(qFunc, c(list(p = up), dparams))
			}

			# tail-sensitive confidence bands
			if (bandType == "equal_local_levels") {
				
				global.bounds <- do.call(get_bounds_two_sided,
				                         c(list(alpha = 1 - conf, n = n), bounds.params))
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

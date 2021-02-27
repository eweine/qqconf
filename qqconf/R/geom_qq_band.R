#' @rdname stat_qq_band
#'
#' @param stat statistic to use to calculate confidence bands. Should be
#'   `qq_band`.
#'
#' @export
geom_qq_band <- function(
	mapping = NULL,
	data = NULL,
	stat = "qq_band",
	position = "identity",
	na.rm = TRUE,
	show.legend = NA,
	inherit.aes = TRUE,
	distribution = "norm",
	dparams = list(),
	difference = FALSE,
	method = "ell",
	alpha = .05
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

  conf <- 1 - alpha
  
	if (conf < 0 | conf > 1) {
		stop("Please provide a valid alpha level for the bands: ",
				 "'conf' must be between 0 and 1.",
				 call. = FALSE)
	}

	# vector with common discrete distributions
	discreteDist <- c("binom", "geom", "nbinom", "pois")

	method <- match.arg(method, c("ell", "ks", "pointwise"))

	ggplot2::layer(
		data = data,
		mapping = mapping,
		stat = stat,
		geom = GeomQqBand,
		position = position,
		show.legend = show.legend,
		inherit.aes = inherit.aes,
		params = list(
			na.rm = na.rm,
			distribution = distribution,
			dparams = dparams,
			difference = difference,
			identity = identity,
			method = method,
			alpha = alpha,
			discrete = distribution %in% discreteDist,
			...
		)
	)
}

#' GeomQqBand
#'
#' @keywords internal
#' @usage NULL
#' @export
GeomQqBand <- ggplot2::ggproto(
	`_class` = "GeomQqBand",
	`_inherit` = ggplot2::Geom,

	default_aes = ggplot2::aes(
		width = 0.75,
		linetype = "solid",
		fontsize = 5,
		shape = 19,
		colour = NA,
		size = .1,
		fill = "blue",
		alpha = .8,
		stroke = 0.1,
		linewidth = .1,
		weight = 1,
		x = NULL,
		y = NULL,
		conds = NULL
	),

	required_aes = c("x", "ymin", "ymax"),

	setup_data = function(data, params) {
		data
	},

	draw_group = ggplot2::GeomRibbon$draw_group,

	draw_key = ggplot2::draw_key_polygon
)

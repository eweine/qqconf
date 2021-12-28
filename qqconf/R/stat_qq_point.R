#' Quantile-quantile points
#'
#' Draws quantile-quantile points, with an additional detrend option.
#'
#' @import ggplot2
#' @importFrom MASS fitdistr
#'
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_point
#'
#' @param distribution The quantile function for the specified distribution. Defaults to qnorm.
#' Custom distributions are allowed as long as all parameters are supplied in dparams.
#' @param difference Whether to plot the difference between the observed and
#'   expected values on the vertical axis.
#' @param log10 Whether to plot axes on -log10 scale (e.g. to see small p-values). Can only be used for strictly
#' positive distributions.
#' @param right_tail This parameter is only used if \code{log10} is \code{TRUE}. When \code{TRUE},
#' the x-axis is -log10(1 - Expected Quantile) and the y-axis is -log10(1 - Observed Quantile).
#' When \code{FALSE} (default) the x-axis is -log10(Expected Quantile) and the y-axis is
#' -log10(Observed Quantile). The parameter should be set to \code{TRUE} to make
#' observations in the right tail of the distribution easier to see, and set to false to make the
#' observations in the left tail of the distribution easier to see.
#' @param dparams List of additional parameters for the quantile function of the distribution
#'   (e.g. df=1). Note that if any parameters of the distribution are specified, parameter estimation will not be performed
#'   on the unspecified parameters, and instead they will take on the default values set by the distribution function.
#'   For the uniform distribution, parameter estimation is not performed, and
#'   the default parameters are max = 1 and min = 0.
#'   For other distributions parameters will be estimated if not provided.
#'   For the normal distribution, we estimate the mean as the median and the standard deviation as \eqn{Sn} from the paper by Rousseeuw and Croux 1993
#'   "Alternatives to the Median Absolute Deviation". For all other distributions besides uniform and normal,
#'   the code uses MLE to estimate the parameters. Note that estimation is not implemented for custom distributions, so all
#'   parameters of the distribution must be provided by the user.
#'
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
  distribution = qnorm,
  dparams = list(),
  difference = FALSE,
  log10 = FALSE,
  right_tail = FALSE,
  ...
) {

  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = StatQqPoint,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      distribution = as.character(substitute(distribution)),
      dparams = dparams,
      difference = difference,
      log10 = log10,
      right_tail = right_tail,
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

  optional_aes = c("label"),

  compute_group = function(data,
                           self,
                           scales,
                           distribution,
                           dparams,
                           difference,
                           log10,
                           right_tail) {

    # distributional function
    qFunc <- distribution

    oidx <- order(data$sample)
    smp <- data$sample[oidx]
    n <- length(smp)
    quantiles <- ppoints(n, a = 0)

    # automatically estimate parameters with MLE
    browser()
    if(length(dparams) == 0) {

      dparams <- estimate_dparams(distribution, smp)

    }

    theoretical <- do.call(qFunc, c(list(p = quantiles), dparams))

    # Here, want to do the same operations as qqplotr on our code
    if (log10 && right_tail) {

      theoretical <- -log10(1 - theoretical)
      if (any(smp <= 0)) {

        stop("log10 scaling can only be used with strictly positive distributions. Consider using a pp plot.")

      }
      smp <- -log10(1 - smp)

    }
    else if (log10 == TRUE) {

      theoretical <- -log10(theoretical)
      if (any(smp <= 0)) {

        stop("log10 scaling can only be used with strictly positive distributions. Consider using a pp plot.")

      }
      smp <- -log10(smp)

    }

    if (difference) {
      smp <- smp - theoretical
    }

    out <- data.frame(sample = smp, theoretical = theoretical)

    if (!is.null(data$label)) out$label <- data$label[oidx]
    out
  }
)

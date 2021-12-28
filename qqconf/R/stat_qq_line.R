#' Quantile-quantile lines
#'
#' Draws a quantile-quantile line, with an additional detrend option.
#'
#' @import ggplot2
#' @importFrom MASS fitdistr
#'
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_path
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
#' @export
stat_qq_line <- function(
  mapping = NULL,
  data = NULL,
  geom = "path",
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
    stat = StatQqLine,
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
             difference,
             log10,
             right_tail) {
      
      # distributional function
      qFunc <- distribution
      
      smp <- sort(data$sample)
      n <- length(smp)
      quantiles <- ppoints(n, a = 0)
      c <- .5 
      
      # automatically estimate parameters with MLE
      if(length(dparams) == 0) {

        dparams <- estimate_dparams(distribution, smp)
        
      }
      
      theoretical <- do.call(qFunc, c(list(p = quantiles), dparams))
      
      if (log10 && right_tail) {
        
        theoretical <- -log10(1 - theoretical)
        low_exp_pt <- c * -log10(do.call(distribution, c(list(p=1 - c(1 / max(n * 1.25, n + 2))), dparams))) + (1 - c) * theoretical[1]
        high_exp_pt <- c * -log10(do.call(distribution, c(list(p=c(1 / max(n * 1.25, n + 2))), dparams))) + (1 - c) * theoretical[n]
        if (any(smp <= 0)) {
          
          stop("log10 scaling can only be used with strictly positive distributions. Consider using a pp plot.")
          
        }
        smp <- -log10(1 - smp)
        
      }
      else if (log10 == TRUE) {
        
        theoretical <- -log10(theoretical)
        low_exp_pt <- c * -log10(do.call(distribution, c(list(p=c(1 / max(n * 1.25, n + 2))), dparams))) + (1 - c) * theoretical[1]
        high_exp_pt <- c * -log10(do.call(distribution, c(list(p=1 - c(1 / max(n * 1.25, n + 2))), dparams))) + (1 - c) * theoretical[n]
        if (any(smp <= 0)) {
          
          stop("log10 scaling can only be used with strictly positive distributions. Consider using pp_conf_plot.")
          
        }
        smp <- -log10(smp)
        
      }
      else {
        
        low_exp_pt <- c * do.call(distribution, c(list(p=c(1 / max(n * 1.25, n + 2))), dparams)) + (1 - c) * theoretical[1]
        high_exp_pt <- c * do.call(distribution, c(list(p=1 - c(1 / max(n * 1.25, n + 2))), dparams)) + (1 - c) * theoretical[n]
        
      }
      
      theoretical <- c(low_exp_pt, theoretical, high_exp_pt)
      
      if (difference) {
        out <- data.frame(xline = c(min(theoretical), max(theoretical)))
        out$yline <- 0
      } else {
        out <- data.frame(xline = c(min(theoretical), max(theoretical)), yline = c(min(theoretical), max(theoretical)))
      }
      
      out$size <- .8
      
      out
    }
  }
)

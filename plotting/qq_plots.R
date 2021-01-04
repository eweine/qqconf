##' QQ plot with global and pointwise confidence bounds
##'
##' Create a qq-plot with with a shaded global confidence interval and,
##' optionally, lines for point-wise bounds. The observed values are plotted
##' against their expected values had they come from the specified distribution.
##'
##' For $N$ independent tests that follow a specified distribution, the global
##' confidence interval of size $C$ is the region such that the probability of
##' all the points being in the region is $C$. The pointwise bounds are determined
##' by taking the marginal probability that an individual point is within the bounds
##' is $C$. If 'difference' is set to TRUE, the vertical axis plots the value
##' Observed quantile minus Expected quantile. Set pw.lty = 0 to suppress plotting
##' of the pointwise bounds.
##' 
##' @param obs The observed data.
##' @param dist The quantile distribution function, qunif by default.
##' @param conf.int The size of the confidence interval for the expected quantiles.
##' @param difference Whether to plot the difference between the observed and
##'   expected values on the vertical axis.
##' @param log10 Whether to plot axes on -log10 scale (e.g. to see small p-values).
##' @param shade.col What color to use for the global confidence interval.
##' @param xlab,ylab Labels for the axes.
##' @param add Whether to add points to an existing plot.
##' @param dist.params List of additional parameters for the quantile distribution
##'   function (e.g. df=1).
##' @param bounds.params List of optional parameters for get_bounds_two_sided
##'   (i.e. tol, max_it, method).
##' @param pw.lty Line type for the pointwise error bounds. Set to zero for no line.
##' @param pw.col Color for the pointwise bounds line.
##' @param ... Additional parameters for the plot.
##'
##' @examples
##' x <- rchisq(1000, 1)
##' qq_conf_plot(x, qchisq, dist.params=list(df=1)) # Plots x against a 1-df chisquare
##'
##' y <- runif(893)
##' qq_conf_plot(y, difference = TRUE, log10 = TRUE, bounds.params = list(method = "search"), pch=3)
qq_conf_plot <- function(obs,
                         dist = qunif,
                         conf.int = 0.95,
                         difference = FALSE,
                         log10 = FALSE,
                         shade.col = 'gray',
                         xlab = 'Expected quantiles',
                         ylab = 'Observed quantiles',
                         add = FALSE,
                         dist.params = NULL,
                         bounds.params = NULL,
                         ## main = NULL,
                         ## ylim = NULL,
                         pw.lty = 3,
                         pw.col = 'black',
                         ...) {
    if(difference && ylab == 'Observed quantiles') {
        ylab <- 'Observed minus expected quantiles'
    }
    samp.size <- length(obs)
    alpha <- signif(1 - conf.int, 4) # Use signif to work around float test in
                                        # get_bound_two_sided
    conf <- c(alpha / 2, conf.int + alpha / 2)
    ## The observed and expected quantiles. Expected quantiles are based on the specified
    ## distribution.
    obs.pts <- sort(obs)
    exp.pts <- do.call(dist, c(list(p=ppoints(samp.size, a=0)), dist.params))
    if (log10 == TRUE) {
        exp.pts <- -log10(exp.pts)
        obs.pts <- -log10(obs.pts)
    }
    if (difference) {
        y.pts <- obs.pts - exp.pts
    } else {
        y.pts <- obs.pts
    }

    ## When not adding points to a qq-plot compute pointwise and global confidence bounds.
    if (!add) {
        left <- exp.pts[1]
        right <- exp.pts[samp.size]
        bottom <- min(y.pts) #obs.pts[1]
        top <- max(y.pts) #obs.pts[samp.size]
        plot(c(left, right), c(bottom, top), type='n', xlab=xlab, ylab=ylab, ...)
        ## plot(c(left, right), c(bottom, top), type = 'n', ...)
        pointwise.low <- do.call(dist,
                                 c(list(p=qbeta(conf[1], 1:samp.size, samp.size:1)), dist.params))
        pointwise.high <- do.call(dist,
                                  c(list(p=qbeta(conf[2], 1:samp.size, samp.size:1)), dist.params))

        ## global.bounds <- get_bounds_two_sided(alpha, samp.size)
        global.bounds <- do.call(get_bounds_two_sided,
                                 c(list(alpha = alpha, n = samp.size), bounds.params))
        global.low <- do.call(dist, c(list(p = global.bounds$lower_bound), dist.params))
        global.high <- do.call(dist, c(list(p = global.bounds$upper_bound), dist.params))
        
        if (log10 == TRUE) {
            pointwise.low <- -log10(pointwise.low)
            pointwise.high <- -log10(pointwise.high)
            global.low <- -log10(global.low)
            global.high <- -log10(global.high)
        }
        if (difference) {
            polygon(c(exp.pts, exp.pts[samp.size:1]),
                    c(global.low - exp.pts, global.high[samp.size:1] - exp.pts[samp.size:1]),
                    border=NA, col=shade.col)
            lines(exp.pts, pointwise.low - exp.pts, lty = pw.lty, col = pw.col, ...)
            lines(exp.pts, pointwise.high - exp.pts, lty = pw.lty, col = pw.col, ...)
        } else {
            polygon(c(exp.pts, exp.pts[samp.size:1]),
                    c(global.low, global.high[samp.size:1]),
                    border=NA, col=shade.col)
            lines(exp.pts, pointwise.low, lty = pw.lty, col = pw.col, ...)
            lines(exp.pts, pointwise.high, lty = pw.lty, col = pw.col, ...)
        }
    }
    points(exp.pts[1:samp.size], y.pts, ...)
    if (difference) {
        abline(h = 0, ...)
    } else {
        abline(0, 1, ...)
    }

    
}


##' QQ plot with global and pointwise confidence bounds
##' 
##' qq-plot with shading of a confidence interval. By default the
##' confidence interval is 95%. Note that the confidence interval shading is
##' only a rough visual aid and not a real confidence interval. It does not,
##' for instance, take into account the correlation of the points in a qq-plot,
##' and it computes the confidence interval for a uniform distribution then
##' transforms those values to the quantiles of the desired distribution, but
##' I don't know that this is the correct confidence interval for the quantiles
##' of the desired distribution.
##' @param obs The observed data
##' @param dist The quantile distribution function, qchisq by default
##' @param conf.int The size of the confidence interval for the expected quantiles
##' @param log10 Whether to plot axes on -log10 scale (e.g. to see small p-values)
##' @param shadecol What color to use for the confidence interval
##' @param add whether to add points to an existing plot
##' @param dist.params: List of additional parameters for the quantile distribution function (e.g. df=1)
##' @param bounds.params List of optional parameters for get_bounds_two_sided (i.e. tol, max_it, method)
##' @param lty Line type for the pointwise error bounds
##' @param ...: additional parameters for the plot
##' Example: qq.conf(x, qchisq, dist.params=list(df=1)) # Plots x against a 1-df chisquare
##' Example: qq.conf(x, qnorm, pch=3) # Plots x against a standard normal
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
        abline(0,1, ...)
    }

    
}

qq_diff_plot <- function(obs, dist=qunif, conf=c(0.025, 0.975), log10=FALSE, shade.col='gray',
                         xlab='Expected quantiles',
                         ylab='Observed - Expected quantiles',
                         add=FALSE,
                     dist.params=NULL,
                     ##pch=20,
                     main=NULL,
                     ylim=NULL,
                     ...) {
    ## qq.conf2 is the same as qq.conf but subtracts the expected from the observed values.
    ## qq-plot with shading of a confidence interval. By default the
    ## confidence interval is 95%. Note that the confidence interval shading is
    ## only a rough visual aid and not a real confidence interval. It does not,
    ## for instance, take into account the correlation of the points in a qq-plot,
    ## and it computes the confidence interval for a uniform distribution then
    ## transforms those values to the quantiles of the desired distribution, but
    ## I don't know that this is the correct confidence interval for the quantiles
    ## of the desired distribution.
    ## obs: The observed data
    ## dist: The quantile distribution function, qchisq by default
    ## conf: The upper and lower bound confidence region to shade
    ## log10: Whether to plot axes on -log10 scale (e.g. to see small p-values)
    ## shade.col: What color to use for the confidence
    ## add: whether to add points to an existing plot
    ## dist.params: List of additional parameters for the quantile distribution function (e.g. df=1)
    ## ...: additional parameters for the plot
    ## Example: qq.conf(x, qchisq, dist.params=list(df=1)) # Plots x against a 1-df chisquare
    ## Example: qq.conf(x, qnorm, pch=3) # Plots x against a standard normal

    samp.size <- length(obs)
    obs.pts <- sort(obs)
    exp.pts <- do.call(dist, c(list(p=ppoints(samp.size, a=0)), dist.params))
    if (log10 == TRUE) {
        exp.pts <- -log10(exp.pts)
        obs.pts <- -log10(obs.pts)
    }
    y.pts <- obs.pts - exp.pts
    if (!add) {
        left <- exp.pts[1]
        right <- exp.pts[samp.size]
        bottom <- min(y.pts)
        top <- max(y.pts)
        plot(c(left, right), c(bottom, top), type='n', xlab=xlab, ylab=ylab, main=main, ylim=ylim)
        e.low <- do.call(dist, c(list(p=qbeta(conf[1], 1:(samp.size), (samp.size):1)), dist.params))
        e.high <- do.call(dist, c(list(p=qbeta(conf[2], 1:(samp.size), (samp.size):1)), dist.params))
        if (log10 == TRUE) {
            e.low <- -log10(e.low)
            e.high <- -log10(e.high)
        }
        polygon(c(exp.pts, exp.pts[(samp.size):1]),
                c(e.low - exp.pts, e.high[(samp.size):1] - exp.pts[samp.size:1]), border=NA, col=shade.col)
    }
    points(exp.pts[1:samp.size], y.pts, ...)
    abline(h=0)
}

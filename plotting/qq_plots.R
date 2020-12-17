##' QQ plot with global confidence bounds
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
##' @param conf The upper and lower bound confidence region to shade
##' @param log10 Whether to plot axes on -log10 scale (e.g. to see small p-values)
##' @param shadecol What color to use for the confidence
##' @param add whether to add points to an existing plot
##' @param distparams: List of additional parameters for the quantile distribution function (e.g. df=1)
##' @param ...: additional parameters for the plot
##' Example: qq.conf(x, qchisq, dist.params=list(df=1)) # Plots x against a 1-df chisquare
##' Example: qq.conf(x, qnorm, pch=3) # Plots x against a standard normal
qq.conf <- function(obs, dist=qunif, conf=c(0.025, 0.975), log10=FALSE, shade.col='gray',
                    xlab='Expected quantiles',
                    ylab='Observed quantiles',
                    add=FALSE,
                    dist.params=NULL,
                    ##pch=20,
                    main=NULL,
                    ylim=NULL,
                    ...) {

    samp.size <- length(obs)
    obs.pts <- sort(obs)
    exp.pts <- do.call(dist, c(list(p=ppoints(samp.size, a=0)), dist.params))
    if (log10 == TRUE) {
        exp.pts <- -log10(exp.pts)
        obs.pts <- -log10(obs.pts)
    }
    if (!add) {
        left <- exp.pts[1]
        right <- exp.pts[samp.size]
        bottom <- obs.pts[1]
        top <- obs.pts[samp.size]
        plot(c(left, right), c(bottom, top), type='n', xlab=xlab, ylab=ylab, main=main, ylim=ylim)
        e.low <- do.call(dist, c(list(p=qbeta(conf[1], 1:(samp.size), (samp.size):1)), dist.params))
        e.high <- do.call(dist, c(list(p=qbeta(conf[2], 1:(samp.size), (samp.size):1)), dist.params))
        if (log10 == TRUE) {
            e.low <- -log10(e.low)
            e.high <- -log10(e.high)
        }
        polygon(c(exp.pts, exp.pts[(samp.size):1]), c(e.low, e.high[(samp.size):1]), border=NA, col=shade.col)
    }
    points(exp.pts[1:samp.size], obs.pts, ...)
    abline(0,1)
}

qq.conf2 <- function(obs, dist=qunif, conf=c(0.025, 0.975), log10=FALSE, shade.col='gray',
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

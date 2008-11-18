`plot.MOStest` <-
    function(x, which = c(1,2,3,6),  ...)
{

    show <- rep(FALSE, 8)
    show[which] <- TRUE
    if (show[1]) {
        X <- x$mod$model$x
        Y <- x$mod$y
        xx <- seq(min(X), max(X), len=101)
        pre <- predict(x$mod, newdata=list(x = xx), se=TRUE)
        g <- x$mod$family$linkinv
        fv <- g(pre$fit)
        hi <- g(pre$fit + 2*pre$se)
        lo <- g(pre$fit - 2*pre$se)
        plot(X, Y, ...)
        matlines(xx, cbind(fv, hi, lo), lty=c(1, 2, 2), lwd=c(2, 1, 1), col=1, ...)
    }
    if (show[2]) {
        require(ellipse) || stop("requires package 'ellipse'")
        ci <- ellipse(x$mod, which=c(2,3))
        plot(ci, type="l", lwd=2, xlim=range(ci[,1],0), ylim=range(ci[,2],0), ...)
        abline(h=0, lty=2, ...)
        par <- x$hump[c("min", "max")]
        par[par==0] <- sqrt(.Machine$double.eps)
        abline(0, -1/2/par[1], ...)
        abline(0, -1/2/par[2], ...)
        k <- coef(summary(x$mod))[2:3, 1:2]
        segments(k[1,1] - k[1,2]*2, k[2,1], k[1,1]+k[1,2]*2, k[2,1], lty=3)
        segments(k[1,1], k[2,1]-k[2,2]*2, k[1,1], k[2,1]+k[2,2]*2, lty=3)
    }
    if (any(show[-c(1,2)])) {
        still <- which(show[-c(1,2)])
        plot(x$mod, which = still, ...)
    }
    invisible()
}

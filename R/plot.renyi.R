`plot.renyi` <-
    function(x, ...)
{
    if (inherits(x, "data.frame")) {
        plt <- factor(rep(rownames(x), ncol(x)), levels=rownames(x))
        alp <- factor(rep(colnames(x), each=nrow(x)), levels=colnames(x))
        div  <- as.vector(as.matrix(x))
        df <- data.frame(diversity=div, plot=plt, alpha=alp)
        lo <- tapply(div, alp, min)
        hi <- tapply(div, alp, max)
        med <- tapply(div, alp, median)
    } else {
        df <- data.frame(diversity = x, alpha = factor(names(x), levels=names(x)), plot = "plot")
        lo <- hi <- med <- NA
    }
    cl <- trellis.par.get("superpose.line")$col
    bwplot(diversity ~ alpha | plot, data=df,  
           panel = function(x, y, ...) {
               panel.lines(x, lo, lty=2, col=cl[3])
               panel.lines(x, med, lty=2, col=cl[2])
               panel.lines(x, hi, lty=2, col=cl[3])
               panel.xyplot(x, y,  ...)
           },
           ...)
}


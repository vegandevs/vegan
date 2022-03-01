`lines.procrustes` <-
    function(x, type=c("segments", "arrows"),  choices=c(1,2),
             truemean = FALSE, ...)
{
    type <- match.arg(type)
    X <- x$X[,choices, drop=FALSE]
    Y <- x$Yrot[, choices, drop=FALSE]
    if (truemean) {
        X <- sweep(X, 2, x$xmean[choices], "+")
        Y <- sweep(Y, 2, x$xmean[choices], "+")
    }
    if (type == "segments")
        ordiArgAbsorber(X[,1], X[,2], Y[,1], Y[,2], FUN = segments, ...)
    else
        ordiArgAbsorber(X[,1], X[,2], Y[,1], Y[,2], FUN = arrows, ...)
    invisible()
}

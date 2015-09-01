"lines.procrustes" <-
    function(x, type=c("segments", "arrows"),  choices=c(1,2), ...)
{
    type <- match.arg(type)
    X <- x$X[,choices, drop=FALSE]
    Y <- x$Yrot[, choices, drop=FALSE]
    if (type == "segments")
        segments(X[,1], X[,2], Y[,1], Y[,2], ...)
    else
        arrows(X[,1], X[,2], Y[,1], Y[,2], ...)
    invisible()
}

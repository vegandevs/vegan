`varpart` <-
    function (Y, X, ..., data, transfo, scale = FALSE) 
{
    if (missing(data)) 
        data <- parent.frame()
    X <- list(X, ...)
    if ((length(X) < 2 || length(X) > 4)) 
        stop("needs 2 to 4 explanatory tables")
    if (!missing(transfo)) {
        Y <- decostand(Y, transfo)
        transfo <- attr(Y, "decostand")
    }
    if (!missing(transfo) && (is.null(dim(Y)) || ncol(Y) == 1)) 
        warning("Transformations probably are meaningless to a single variable")
    if (scale && !missing(transfo)) 
        warning("Y should not be both transformed and scaled (standardized)")
    Y <- scale(Y, center = TRUE, scale = scale)
    Sets <- list()
    for (i in 1:length(X)) {
        if (inherits(X[[i]], "formula")) {
            mf <- model.frame(X[[i]], data, na.action = na.fail,
                              drop.unused.levels = TRUE)
            trms <- attr(mf, "terms")
            Sets[[i]] <- model.matrix(trms, mf)
            if (any(colnames(Sets[[i]]) == "(Intercept)")) {
                xint <- which(colnames(Sets[[i]]) == "(Intercept)")
                Sets[[i]] <- (Sets[[i]])[, -xint, drop = FALSE]
            }
        }
        else Sets[[i]] <- as.matrix(X[[i]])
        Sets[[i]] <- scale(Sets[[i]], center = TRUE, scale = TRUE)
    }
    out <- list()
    out$part <- switch(length(Sets), NULL,
                       varpart2(Y, Sets[[1]], Sets[[2]]),
                       varpart3(Y, Sets[[1]], Sets[[2]], Sets[[3]]), 
                       varpart4(Y, Sets[[1]], Sets[[2]], Sets[[3]], Sets[[4]]))
    out$scale <- scale
    if (!missing(transfo)) 
        out$transfo <- transfo
    out$call <- match.call()
    mx <- rep(" ", length(X))
    for (i in 1:length(X)) mx[i] <- deparse(out$call[[i+2]], width.cutoff = 500)
    out$tables <- mx
    class(out) <- c("varpart", class(out))
    out
}

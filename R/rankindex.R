"rankindex" <-
function (grad, veg, indices = c("euc", "man", "gow", "bra", 
    "kul"), stepacross = FALSE, method = "spearman", ...) 
{
    grad <- as.data.frame(grad)
    if (any(sapply(grad, is.factor))) {
        require(cluster) || stop("factors in 'grad' need package 'cluster'")
        message("'grad' included factors: used cluster:::daisy")
        span <- daisy(grad)
    } else {
        span <- vegdist(grad, "eucl")
    }
    veg <- as.matrix(veg)
    res <- numeric(length(indices))
    ## create names if indices is a list of functions without names
    if (is.list(indices)) {
        nam <- names(indices)
        if (is.null(nam))
            nam <- paste("dis", 1:length(indices), sep="")
    } else nam <- indices
    names(res) <- nam
    ## indices is a list of functions which return dist objects
    if (is.list(indices)) {
        for (i in 1:length(indices)) {
            ## don't accept similarities
            if (indices[[i]](matrix(1, 2, 2)) != 0)
                stop("define dissimilarity and not similarity")
            y <- indices[[i]](veg)
            ## check class of output
            if (!inherits(y, "dist"))
                stop("function in 'indices' must return a 'dist' object")
            if (stepacross) {
                is.na(y) <- no.shared(veg)
                y <- stepacross(y, trace = FALSE, toolong = -1, ...)
            }
            res[i] <- cor(span, y, method = method)
        }
    ## indices is a character vector naming methods in vegdist
    } else {
        for (i in indices) {
            y <- vegdist(veg, i)
            if (stepacross) {
                is.na(y) <- no.shared(veg)
                y <- stepacross(y, trace = FALSE, toolong = -1, ...)
            }
            res[i] <- cor(span, y, method = method)
        }
    }
    res
}

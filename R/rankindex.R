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
    names(res) <- indices
    for (i in indices) {
        y <- vegdist(veg, i)
        if (stepacross) {
            is.na(y) <- no.shared(veg)
            y <- stepacross(y, trace = FALSE, toolong = -1, ...)
        }
        res[i] <- cor(span, y, method = method)
    }
    res
}


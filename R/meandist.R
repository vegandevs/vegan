`meandist` <-
    function(dist, grouping, ...)
{
    ## merge levels so that lower is always first (filling lower triangle)
    mergenames <- function(X, Y, ...) {
        xy <- cbind(X, Y)
        xy <- apply(xy, 1, sort)
        apply(xy, 2, paste, collapse = " ")
    }
    grouping <- factor(grouping, exclude = NULL)
    cl <- outer(grouping, grouping, mergenames)
    cl <- cl[lower.tri(cl)]
    ## Cannot have within-group dissimilarity for group size 1
    n <- table(grouping)
    take <- matrix(TRUE, nlevels(grouping), nlevels(grouping))
    diag(take) <- n > 1
    take[upper.tri(take)] <- FALSE
    ## Get output matrix
    out <- matrix(NA, nlevels(grouping), nlevels(grouping))
    out[take] <- tapply(dist, cl, mean)
    out[upper.tri(out)] <- t(out)[upper.tri(out)]
    rownames(out) <- colnames(out) <- levels(grouping)
    class(out) <- c("meandist", "matrix")
    attr(out, "n") <- table(grouping)
    out
}


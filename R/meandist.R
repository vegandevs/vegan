`meandist` <-
    function(dist, grouping, ...)
{
    if (!inherits(dist, "dist"))
        stop("'dist' must be dissimilarity object inheriting from", dQuote(dist))
    ## check that 'dist' are dissimilarities (non-negative)
    if (any(dist < -sqrt(.Machine$double.eps)))
        warning("some dissimilarities are negative: is this intentional?")
    grouping <- factor(grouping, exclude = NULL)
    ## grouping for rows and columns
    grow <- grouping[as.dist(row(as.matrix(dist)))]
    gcol <- grouping[as.dist(col(as.matrix(dist)))]
    ## The row index must be "smaller" of the factor levels so that
    ## all means are in the lower triangle, and upper is NA
    first <- as.numeric(grow) >= as.numeric(gcol)
    cl1 <- ifelse(first, grow, gcol)
    cl2 <- ifelse(!first, grow, gcol)
    ## Cannot have within-group dissimilarity for group size 1
    n <- table(grouping)
    take <- matrix(TRUE, nlevels(grouping), nlevels(grouping))
    diag(take) <- n > 1
    take[upper.tri(take)] <- FALSE
    out <- matrix(NA, nlevels(grouping), nlevels(grouping))
    ## Get output matrix
    tmp <- tapply(dist, list(cl1, cl2), mean)
    out[take] <- tmp[!is.na(tmp)]
    out[upper.tri(out)] <- t(out)[upper.tri(out)]
    rownames(out) <- colnames(out) <- levels(grouping)
    class(out) <- c("meandist", "matrix")
    attr(out, "n") <- table(grouping)
    out
}


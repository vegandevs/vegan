### Rotates metaMDS result so that axis one is parallel to vector 'x'
`metaMDSrotate` <-
    function(object, vec, na.rm = FALSE, ...) 
{
    if (!inherits(object, "metaMDS"))
        stop(gettextf("function works only with 'metaMDS' results"))
    x <- object$points
    sp <- object$species
    N <- NCOL(x)
    if (N < 2)
        stop(gettextf("needs at least 2 dimensions"))
    vec <- drop(vec)
    if (length(dim(vec)) > 1)
        stop(gettextf("function works only with univariate 'vec'"))
    if (!is.numeric(vec))
        stop(gettextf("'vec' must be numeric"))
    ## scores must be orthogonal for the next loop to work
    if (N > 2) {
        pc <- prcomp(x)
        x <- pc$x
        if (!all(is.na(sp)))
            sp <- sp %*% pc$rotation
    }
    ## vectorfit finds the direction cosine. We rotate first axis to
    ## 'vec' which means that we make other axes orthogonal to 'vec'
    ## one by one
    if (na.rm)
        keep <- !is.na(vec)
    else
        keep <- !logical(length(vec))
    for (k in 2:N) {
        rot <- vectorfit(x[keep, c(1,k)], vec[keep], permutations=0)$arrows
        rot <- drop(rot)
        ## counterclockwise rotation matrix:
        ## [cos theta   -sin theta]
        ## [sin theta    cos theta]
        rot <- rbind(rot, rev(rot))
        rot[1,2] <- -rot[1,2]
        ## Rotation of points and species scores
        x[, c(1,k)] <- x[, c(1,k)] %*% rot
        if (!all(is.na(sp)))
            sp[, c(1,k)] <- sp[, c(1,k)] %*% rot
    }
    ## Rotate 2..N axes to PC
    if (N > 2 && attr(object$points, "pc")) {
        pc <- prcomp(x[,-1])
        x[,-1] <- pc$x
        if (!all(is.na(sp)))
            sp[,-1] <- sp[,-1] %*% pc$rotation
    }
    ## '[] <-' retains attributes
    object$points[] <- x
    object$species[] <- sp
    attr(object$points, "pc") <- FALSE
    object
}


### Rotates metaMDS result so that axis one is parallel to vector 'x'
`metaMDSrotate` <-
    function(object, vec, choices = 1:2, ...) 
{
    if (!inherits(object, "metaMDS"))
        stop("function works only with 'metaMDS' results")
    if (length(choices) != 2)
        stop("function can be only used with 2dim plots")
    vec <- drop(vec)
    if (length(dim(vec)) > 1)
        stop("function works only with univariate 'x'")
    ## envfit finds the direction cosine
    rot <- envfit(object, vec, choices = choices, ...)$vectors$arrows
    rot <- drop(rot)
    ## rotation matrix [[sin theta, cos theta] [-cos theta, sin theta]]
    rot <- rbind(rot, rev(rot))
    rot[2,1] <- -rot[2,1]
    ## transpose (inverse) of the rotation matrix
    rot <- t(rot)
    ## Rotation of points and species scores
    object$points[] <- object$points %*% rot
    attr(object$points, "pc") <- FALSE
    if (!is.null(object$species))
        object$species[] <- object$species %*% rot
    object
}


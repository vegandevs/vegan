### Rotates metaMDS or monoMDS result so that axis one is parallel to
### vector 'x'.
`MDSrotate` <-
    function(object, vec, na.rm = FALSE, ...)
{
    workswith <- c("metaMDS", "monoMDS", "GO")
    if (!inherits(object, workswith))
        stop(gettextf("function works only with the results of: %s",
                      paste(workswith, collapse = ", ")))
    x <- object$points
    if (is.null(object$species))
        sp <- NA
    else
        sp <- object$species
    N <- NCOL(x)
    if (N < 2)
        stop(gettextf("needs at least two dimensions"))
    ## check if vec is a factor and then use lda to find a matrix that
    ## separates optimally factor levels
    if (is.factor(vec) || is.character(vec)) {
        da <- lda(x, vec)
        vec <- predict(da, dimen = N - 1)$x
        message(sprintf(ngettext(NCOL(vec),
                         "factor replaced with discriminant axis",
                         "factor replaced with %d discriminant axes",
                                 ), NCOL(vec)))
        if (NCOL(vec) > 1)
            message(gettextf("proportional traces: %.3f",
                             da$svd[1:NCOL(vec)]^2/sum(da$svd^2)))
    }
    vec <- as.matrix(vec)
    NV <- NCOL(vec)
    if (NV >= N)
        stop(gettextf("you can have max %d vectors, but you had %d",
             N-1, NV))
    if (!is.numeric(vec))
        stop(gettextf("'vec' must be numeric"))
    ## vectorfit finds the direction cosine. We rotate first axis to
    ## 'vec' which means that we make other axes orthogonal to 'vec'
    ## one by one
    if (na.rm)
        keep <- complete.cases(vec)
    else
        keep <- !logical(NROW(vec))
    ## Rotation loop
    for(v in seq_len(NV)) {
        for (k in (v+1):N) {
            arrs <- vectorfit(x[keep,], vec[keep,v], permutations = 0)$arrows
            rot <- arrs[c(v,k)]/sqrt(sum(arrs[c(v,k)]^2))
            rot <- drop(rot)
            ## counterclockwise rotation matrix:
            ## [cos theta   -sin theta]
            ## [sin theta    cos theta]
            rot <- rbind(rot, rev(rot))
            rot[1,2] <- -rot[1,2]
            ## Rotation of points and species scores
            x[, c(v,k)] <- x[, c(v,k)] %*% rot
            if (!all(is.na(sp)))
                sp[, c(v,k)] <- sp[, c(v,k)] %*% rot
        }
    }
    ## Two or more free axes are (optionally) rotated to PCs
    if (N - NV > 1 && attr(object$points, "pc")) {
        pc <- prcomp(x[,-seq_len(NV)])
        x[,-seq_len(NV)] <- pc$x
        if (!all(is.na(sp)))
            sp[,-seq_len(NV)] <- sp[,-seq_len(NV)] %*% pc$rotation
    }
    ## '[] <-' retains attributes
    object$points[] <- x
    object$species[] <- sp
    attr(object$points, "pc") <- FALSE
    object
}


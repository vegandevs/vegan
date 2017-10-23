`calibrate.cca` <-
    function(object, newdata, rank = "full", ...)
{
    ## inversion solve(b) requires a square matrix, and we should
    ## append imaginary dims to get those in dbrda with negative
    ## constrained eigenvalues. Work is need to to verify this can be
    ## done, and therefore we just disable calibrate with negative
    ## eigenvalues in constraints.
    if (inherits(object, "dbrda") && object$CCA$poseig < object$CCA$qrank)
        stop("cannot be used with 'dbrda' with imaginary constrained dimensions")
    if (!is.null(object$pCCA))
        stop("does not work with conditioned (partial) models")
    if (is.null(object$CCA) || object$CCA$rank == 0)
        stop("needs constrained model")
    if (object$CCA$rank < object$CCA$qrank)
        stop("rank of constraints is higher than rank of dependent data")
    if (rank != "full")
        rank <- min(rank, object$CCA$rank)
    else
        rank <- object$CCA$rank
    if (missing(newdata))
        wa <- object$CCA$wa
    else
        wa <- predict(object, type="wa", newdata=newdata)
    qrank <- object$CCA$qrank
    b <- (coef(object))[object$CCA$QR$pivot[1:qrank], , drop=FALSE]
    b <- solve(b)
    pred <- wa[ , 1:rank, drop=FALSE]  %*% b[1:rank, , drop =FALSE]
    envcen <- object$CCA$envcentre[object$CCA$QR$pivot]
    envcen <- envcen[1:object$CCA$qrank]
    sweep(pred, 2, envcen, "+")
}


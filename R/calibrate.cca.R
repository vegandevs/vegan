`calibrate.cca` <-
    function(object, newdata, rank = "full", ...)
{
    if (!is.null(object$pCCA))
        stop("does not work with conditioned (partial) models")
    if (is.null(object$CCA))
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
    pred <- wa[ , 1:rank, drop=FALSE]  %*% b[1:qrank, , drop =FALSE]
    envcen <- object$CCA$envcentre[object$CCA$QR$pivot]
    envcen <- envcen[1:object$CCA$qrank]
    pred <- sweep(pred, 2, envcen, "+")
    pred
}


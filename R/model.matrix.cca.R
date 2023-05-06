`model.matrix.cca` <-
    function(object, ...)
{
    X <- Z <- NULL
    w <- 1/sqrt(object$rowsum)
    if (!is.null(object$pCCA))
        Z <- w * qr.X(object$pCCA$QR, ncol = length(object$pCCA$QR$pivot))
    if (!is.null(object$CCA)) {
        X <- qr.X(object$CCA$QR, ncol = length(object$CCA$QR$pivot))
        ## First columns come from Z
        if (!is.null(Z))
            X <- X[, -seq_len(ncol(Z)), drop = FALSE]
        X <- w * X
    }
    m <- list()
    if (!is.null(Z))
        m$Conditions <- Z
    if (!is.null(X))
        m$Constraints <- X
    if (length(m) == 1)
        m <- m[[1]]
    m
}

`model.matrix.rda` <-
    function(object, ...)
{
    X <- Z <- NULL
    if (!is.null(object$pCCA))
        Z <- qr.X(object$pCCA$QR, ncol = length(object$pCCA$QR$pivot))
    if (!is.null(object$CCA)) {
        X <- qr.X(object$CCA$QR, ncol = length(object$CCA$QR$pivot))
        if (!is.null(Z))
            X <- X[, -seq_len(ncol(Z)), drop=FALSE]
    }
    m <- list()
    if (!is.null(Z))
        m$Conditions <- Z
    if (!is.null(X))
        m$Constraints <- X
    if (length(m) == 1)
        m <- m[[1]]
    m
}

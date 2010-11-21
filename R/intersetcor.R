intersetcor <- function(object) 
{
    if (!inherits(object, "cca"))
        stop("can be used only with objects inheriting from 'cca'")
    wa <- object$CCA$wa
    if (!inherits(object, "rda")) {   # is CCA
        w <- object$rowsum
        wa <- sweep(object$CCA$wa, 1, sqrt(w), "*")
    }
    X <- qr.X(object$CCA$QR)
    ## current R (2.12.0) uses wrong column names in pivoted qr.X()
    if (getRversion() <= "2.12.0")
        colnames(X)[object$CCA$QR$pivot] <- colnames(X)
    cor(X, wa)
}

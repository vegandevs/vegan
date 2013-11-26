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
    cor(X, wa)
}

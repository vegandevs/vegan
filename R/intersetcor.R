intersetcor <- function(object) 
{
    if (!inherits(object, "cca"))
        stop("can be used only with objects inheriting from 'cca'")
    w <- weights(object)
    wa <- sweep(object$CCA$wa, 1, sqrt(w), "*")
    cor(qr.X(object$CCA$QR), wa)
}

`fitted.capscale` <-
    function(object, model = c("CCA", "CA", "pCCA", "Imaginary"),
             type = c("response", "working"), ...)
{
    model <- match.arg(model)
    type <- match.arg(type)
    ## Return scaled eigenvalues
    U <- switch(model,
                CCA = object$CCA$u %*% diag(sqrt(object$CCA$eig)),
                CA = object$CA$u %*% diag(sqrt(object$CA$eig)),
                Imaginary = object$CA$imaginary.u.eig,
                pCCA = object$pCCA$Fit/object$adjust)
    ## Distances or working scores U
    if (type == "response") {
        U <- dist(U)
        ## remove additive constant (if add = TRUE)
        if (!is.null(object$ac))
            U <- U - object$ac
    }
    U
}

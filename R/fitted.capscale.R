fitted.capscale <-
    function(object, model = c("CCA", "CA", "Imaginary"),
             type = c("response", "working"), ...)
{
    model <- match.arg(model)
    type <- match.arg(type)
    ## Return scaled eigenvalues
    U <- switch(model,
                CCA = object$CCA$u.eig,
                CA = object$CA$u.eig,
                Imaginary = object$CA$imaginary.u.eig)
    ## Distances or working scores U
    if (type == "response")
        dist(U)
    else
        U
}

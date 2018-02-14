`fitted.capscale` <-
    function(object, model = c("CCA", "CA", "pCCA", "Imaginary"),
             type = c("response", "working"), ...)
{
    model <- match.arg(model)
    if (is.null(object[[model]]) && model != "Imaginary")
        stop(gettextf("component '%s' does not exist", model))
    type <- match.arg(type)
    ## Return scaled eigenvalues
    U <- switch(model,
                CCA = ordiYbar(object, "CCA"),
                CA = ordiYbar(object, "CA"),
                Imaginary = object$CA$imaginary.u.eig,
                pCCA = ordiYbar(object, "pCCA"))
    if (is.null(U))
        stop(gettextf("component '%s' does not exist", model))
    ## Distances or working scores U
    if (type == "response") {
        U <- dist(U)
        ## undo sqrt.dist -- sqrt.dist was applied first in capscale,
        ## so it must be last here
        if (object$sqrt.dist)
            U <- U^2
    }
    U
}

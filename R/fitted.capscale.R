`fitted.capscale` <-
    function(object, model = c("CCA", "CA", "pCCA", "Imaginary"),
             type = c("response", "working"), ...)
{
    model <- match.arg(model)
    if (is.null(object[[model]]) && model != "Imaginary")
        stop("component ", model, " does not exist")
    type <- match.arg(type)
    ## Return scaled eigenvalues
    U <- switch(model,
                CCA = ordiYbar(object, "CCA"),
                CA = ordiYbar(object, "CA"),
                Imaginary = object$CA$imaginary.u.eig,
                pCCA = ordiYbar(object, "pCCA"))
    if (is.null(U))
        stop("component ", model, " does not exist")
    ## Distances or working scores U
    if (type == "response") {
        U <- dist(U)
        ## remove additive constant (if add = TRUE)
        if (!is.null(object$ac)) {
            if (object$add == "lingoes")
                U <- sqrt(U^2 - 2 * object$ac)
            else if (object$add == "cailliez")
                U <- U - object$ac
            else
                stop("unknown Euclidifying adjustment")
        }
        ## undo sqrt.dist -- sqrt.dist was applied first in capscale,
        ## so it must be last here
        if (object$sqrt.dist)
            U <- U^2
    }
    U
}

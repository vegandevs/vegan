`summary.cca` <- function (object,
                           digits = max(3, getOption("digits") - 3),
                           ...)
{
    summ <- list()
    summ$call <- object$call
    summ$tot.chi <- object$tot.chi
    ## only the Real component for capscale() with negative eigenvalues
    if (!is.null(object$CA$imaginary.chi))
        summ$tot.chi <- summ$tot.chi - object$CA$imaginary.chi
    summ$partial.chi <- object$pCCA$tot.chi
    summ$constr.chi <- object$CCA$tot.chi
    summ$unconst.chi <- object$CA$tot.chi
    ## nested list cont$importance needed to keep vegan pre-2.5-0 compatibility
    summ$cont$importance <- summary(eigenvals(object))
    if (!is.null(object$CCA) && object$CCA$rank > 0)
        summ$concont$importance <- summary(eigenvals(object, model = "constrained"))
    summ$digits <- digits
    summ$inertia <- object$inertia
    summ$method <- object$method
    class(summ) <- "summary.cca"
    summ
}

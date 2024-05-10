`summary.cca` <- function (object,
                           digits = max(3, getOption("digits") - 3),
                           ...)
{
    ## From vegan 2.6-6 summary.cca no longer spits information on
    ## scores, but many packages still use summary() to access
    ## scores. This code tries to handle those functions so that the
    ## package maintainers have time to switch from summary() to
    ## scores() to get scores. This code will be missing in later
    ## releases and then summary() will fail to get scores.
    dots <- match.call(expand.dots = FALSE)$...
    if (isTRUE(any(names(dots) %in%
                   c("scaling", "axes", "display", "correlation", "hill"))))
        warning("summary() to get scores is deprecated: use scores()")
    if (isTRUE("axes" %in% names(dots)))
        choices <- seq_len(dots$axes)
    else
        choices <- seq_len(6)
    summ <- scores(object, choices = choices, ...)
    ## mercy code of handling scores in summary ends.
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

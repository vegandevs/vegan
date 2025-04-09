`drop1.cca` <-
    function(object, scope, test = c("none", "permutation"),
             permutations = how(nperm = 199), ...)
{
    if (inherits(object, "prc"))
        stop("'step'/'drop1' cannot be used for 'prc' objects")
    if (is.null(object$terms))
        stop("ordination model must be fitted using formula")
    test <- match.arg(test)
    # don't show messages about aliased terms
    out <- suppressMessages(NextMethod("drop1", object, test = "none"))
    cl <- class(out)
    if (test == "permutation") {
        rn <- rownames(out)[-1]
        if (missing(scope))
            scope <- rn
        else if (!is.character(scope))
            scope <- drop.scope(scope)
        adds <- anova(object, by = "margin", scope = scope,
                      permutations = permutations, ...)
        out <- cbind(out, rbind(NA, adds[rn,3:4]))
        class(out) <- cl
    }
    out
}

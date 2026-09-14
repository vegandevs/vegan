`anova.prc` <-
    function(object, ...)
{
    ## if user specified 'by', cast prc() to an rda() and call anova
    ## on its result
    extras <- match.call(expand.dots = FALSE)
    if ("by" %in% names(extras$...)) {
        Y <- deparse1(object$call$response)
        X <- deparse1(object$call$treatment)
        Z <- deparse1(object$call$time)
        fla <- paste(Y, "~", X, "*", Z, "+ Condition(", Z, ")")
        fla <- as.formula(fla)
        ## get extras
        m <- match(c("data", "scale", "subset", "na.action"),
                   names(object$call), 0)
        call <- object$call[c(1,m)]
        call$formula <- fla
        call[[1]] <- as.name("rda")
        object <- eval(call, parent.frame())
        anova(object, ...)
    } else {
        NextMethod("anova")
    }
}

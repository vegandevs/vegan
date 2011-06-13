`anova.prc` <-
    function(object, ...)
{
    ## if user specified 'by', cast prc() to an rda() and call anova
    ## on its result
    extras <- match.call(expand.dots = FALSE)
    if ("by" %in% names(extras$...)) {
        Y <- as.character(object$call$response)
        X <- as.character(object$call$treatment)
        Z <- as.character(object$call$time)
        fla <- paste(Y, "~", X, "*", Z, "+ Condition(", Z, ")")
        fla <- as.formula(fla)
        ## get extras
        m <- match(c("data", "scale", "subset", "na.action"),
                   names(object$call), 0)
        call <- object$call[c(1,m)]
        call$formula <- fla
        call[[1]] <- as.name("rda.formula")
        object <- eval(call, parent.frame())
        anova(object, ...)
    } else {
        NextMethod("anova", object, ...)
    }    
}

"prc" <-
    function (response, treatment, time, ...) 
{
    extras <- match.call(expand.dots=FALSE)$...
    if (is.null(extras$data))
        data <- parent.frame()
    else
        data <- eval(extras$data)
    y <- deparse(substitute(response))
    x <- deparse(substitute(treatment))
    z <- deparse(substitute(time))
    fla <- as.formula(paste("~", x, "+", z))
    mf <- model.frame(fla, data)
    if (!all(sapply(mf, is.factor)))
        stop(x, " and ", z, " must be factors")
    if (any(sapply(mf, is.ordered))) 
        stop(x, " or ", z, " cannot be ordered factors")
    fla <- as.formula(paste(y, "~", z, "*", x, "+ Condition(", 
                            z, ")"))
    mod <- rda(fla, ...)
    mod$call <- match.call()
    class(mod) <- c("prc", class(mod))
    mod
}

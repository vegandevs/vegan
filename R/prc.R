`prc`  <-
    function (response, treatment, time, ...)
{
    extras <- match.call(expand.dots = FALSE)$...
    if (is.null(extras$data))
        data <- parent.frame()
    else
        data <- eval(extras$data)
    y <- deparse(substitute(response))
    x <- deparse(substitute(treatment))
    z <- deparse(substitute(time))
    oldcon <- options(contrasts = c("contr.treatment", "contr.poly"))
    on.exit(options(oldcon))
    fla <- as.formula(paste("~", x, "+", z))
    mf <- model.frame(fla, data, na.action = na.pass)
    if (!all(sapply(mf, is.factor)))
        stop(gettextf("%s and %s must be factors", x, z))
    if (any(sapply(mf, is.ordered)))
        stop(gettextf("%s or %s cannot be ordered factors", x, z))
    fla.zx <- as.formula(paste("~", z, ":", x))
    fla.z <- as.formula(paste("~", z))
    # delete first (control) level from the design matrix
    X = model.matrix(fla.zx, mf)[,-c(seq_len(nlevels(time)+1))]
    Z = model.matrix(fla.z, mf)[,-1]
    mod <- rda(response ~ X + Condition(Z), ...)
    mod$terminfo$xlev = list(levels(time), levels(treatment))
    names(mod$terminfo$xlev) = c(paste(z), paste(x))
    mod$call <- match.call()
    class(mod) <- c("prc", class(mod))
    mod
}

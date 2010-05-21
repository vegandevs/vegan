`ordiresids` <-
    function(x, kind = c("residuals", "scale", "qqmath"), residuals = "working",
             type = c("p", "smooth", "g"), formula, ...)
{
    require(lattice) || stop("requires package lattice")
    kind <- match.arg(kind)
    if (!inherits(x, "cca") || is.null(x$CCA))
        stop("function is only available for constrained ordination")
    fit <- fitted(x, type = residuals)
    res <- residuals(x, type = residuals)
    ## remove the effects of row weights in CA
    if (!inherits(x, "rda")) {
        sqr <- sqrt(x$rowsum)
        fit <- sweep(fit, 1, sqr, "*")
        res <- sweep(res, 1, sqr, "*")
    }
    colnam <- rep(colnames(fit), each=nrow(fit))
    rownam <- rep(rownames(fit), ncol(fit))
    df <- data.frame(Fitted = as.vector(fit), Residuals = as.vector(res))
    if (!is.null(rownam))
        df$Sites <- rownam
    if (!is.null(colnam))
        df$Species <- colnam
    if (kind == "residuals") {
        if (missing(formula))
            formula <- as.formula(Residuals ~ Fitted)
        pl <- xyplot(formula, data = df, type = type, ...)
    }
    if (kind == "scale") {
        if (missing(formula))
            formula <- as.formula(sqrt(abs(Residuals)) ~ Fitted)
        pl <- xyplot(formula, data = df, type = type, ...)
    }
    if (kind == "qqmath") {
        if (missing(formula))
            formula <- as.formula(~ Residuals)
        pl <- qqmath(formula, data = df, type = type, ...)
    }
    pl
}

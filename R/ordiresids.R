`ordiresids` <-
    function(x, kind = c("residuals", "scale", "qqmath"),
             residuals = "working",
             type = c("p", "smooth", "g"), formula, ...)
{
    message("ordiresids is deprecated:\nuse directly fitted, residuals, rstandard, rstudent etc.")
    kind <- match.arg(kind)
    if (!inherits(x, "cca") || is.null(x$CCA) || x$CCA$rank == 0)
        stop("function is only available for constrained ordination")
    residuals <- match.arg(residuals, c("working", "response",
                                        "standardized", "studentized"))
    fit <- switch(residuals,
        "working" =,
        "response" = fitted(x, type = residuals),
        "standardized" =,
        "studentized" = sweep(fitted(x, type="working"), 2, sigma(x), "/"))
    res <- switch(residuals,
        "standardized" = rstandard(x),
        "studentized" = rstudent(x),
        residuals(x, type = residuals))
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

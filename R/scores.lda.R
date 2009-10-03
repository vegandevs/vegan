`scores.lda` <-
    function(x, display, ...)
{
    require(MASS) || stop("'lda' objects created in MASS need MASS for 'scores'")
    display <- match.arg(display,
                         c("sites", "species", "scores", "predictors", "x", "coef"),
                         several.ok = TRUE)
    out <- NULL
    if (display %in% c("sites", "scores", "x"))
        out[["scores"]] <- predict(x)$x
    if (display %in% c("species", "predictors", "coef"))
        out[["coefficients"]] <- coef(x)
    if (length(out) == 1)
        out <- out[[1]]
    out
}

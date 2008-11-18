`radlattice` <-
    function(x, BIC = FALSE, ...)
{
    require(lattice) || stop("requires package 'lattice'")
    y <- x$y
    fv <- unlist(fitted(x))
    mods <- names(x$models)
    p <- length(mods)
    n <- length(y)
    Abundance <- rep(y, p)
    Rank <- rep(1:n, p)
    Model <- factor(rep(mods, each=n), levels = mods)
    aic <- AIC(x, BIC = BIC)
    col <- trellis.par.get("superpose.line")$col
    if (length(col) > 1)
        col <- col[2]
    xyplot(Abundance ~  Rank | Model, subscripts = TRUE,
           scales = list(y = list(log = 2)), as.table = TRUE,
           panel = function(x, y, subscripts) {
               panel.xyplot(x, y, ...)
               panel.xyplot(x, log2(fv[subscripts]), type="l", lwd=3,
                            col = col, ...)
               panel.text(max(x), max(y), paste(if (BIC) "BIC" else "AIC", "=",
                          formatC(aic[panel.number()], digits=2, format="f")),
                          pos=2)
           }
           )
}


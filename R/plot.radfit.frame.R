`plot.radfit.frame` <-
    function (x, which, BIC = FALSE, model, legend = TRUE, ...)
{
    modnam <- names(x[[1]]$models)
    if (!missing(model))
        pick <- pmatch(model, modnam, nomatch = FALSE)
    else pick <- FALSE
    modcols <- seq_along(modnam)
    if (missing(which))
        which <- seq_len(length(x))
    npanels <- prod(par("mfrow"))
    if (npanels < length(which))
        warning("you have ", length(x), " models, but graph only has ", npanels,
                " panels")
    for (panel in which) {
        plot(x[[panel]]$y, ...)
        fv <- fitted(x[[panel]])
        if (pick) {
            take <- pick
        } else {
            if (BIC)
                k <- log(length(x[[panel]]$y))
            else
                k <- 2
            take <- which.min(AIC(x[[panel]], k))
        }
        fv <- fv[, take, drop=FALSE]
        lines(fv, col = modcols[take], ...)
        title(main = names(x)[panel])
        if (legend)
            legend("topright", modnam[take], col=modcols[take], lty=1, ...)
    }
}

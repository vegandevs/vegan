### support functions for wcmdscale results: print, scores and plot.

`print.wcmdscale` <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{
    writeLines(strwrap(pasteCall(x$call)))
    cat("\n")
    ## tabulate total inertia and ranks
    totev <- sum(x$eig)
    negax <- x$eig < 0
    if (any(negax)) {
        ranks <- c(NA, sum(!negax), sum(negax))
        negax <- x$eig < 0
        realev <- sum(x$eig[!negax])
        imev <- sum(x$eig[negax])
        evs <- c("Total" = totev, "Real" = realev, "Imaginary" = imev)
    } else {
        ranks <- length(x$eig)
        evs <- c("Total" = totev)
    }
    tbl <- cbind("Inertia" = evs, "Rank" = ranks)
    printCoefmat(tbl, digits = digits, na.print = "")
    if (!is.na(x$ac) && x$ac > 0)
        cat("additive constant ", x$ac, " (method ", x$add, ")\n", sep = "")
    cat("\nResults have", NROW(x$points), "points,", NCOL(x$points), "axes\n")
    ## print eigenvalues, but truncate very long lists
    PRINLIM <- 120
    neig <- length(x$eig)
    cat("\nEigenvalues:\n")
    print(zapsmall(x$eig[1 : min(neig, PRINLIM)], digits = digits, ...))
    if (neig > PRINLIM)
        cat("(Showing", PRINLIM, "of", neig, "eigenvalues)\n")
    wvar <- var(x$weights)
    wlen <- length(x$weights)
    cat("\nWeights:")
    if (wvar < 1e-6)
        cat(" Constant\n")
    else {
        cat("\n")
        print(zapsmall(x$weights[1 : min(wlen, PRINLIM)], digits = digits, ...))
        if (wlen > PRINLIM)
            cat("(Showing", PRINLIM, "of", wlen, "weights)\n")
    }
    cat("\n")
    invisible(x)
}

`scores.wcmdscale` <-
    function(x, choices = NA, display = "sites", tidy = FALSE, ...)
{
    if ("species" %in% names(x))
        display <- match.arg(display, c("sites", "species"), several.ok = TRUE)
    else
        display <- match.arg(display, "sites")
    p <- list()
    for (sco in display) {
        if (sco == "sites")
            p[[sco]] <- x$points
        else
            p[[sco]] <- x[[sco]]
        if (!anyNA(choices)) {
            choices <- choices[choices <= NCOL(p[[sco]])]
            p[[sco]] <- p[[sco]][, choices, drop = FALSE]
        }
        if (is.null(rownames(p[[sco]])))
            rownames(sco) <- as.character(seq_len(nrow(p[[sco]])))
        if (tidy) {
            wts <- if(sco == "sites") weights(x) else NA
            p[[sco]] <- data.frame(p[[sco]], "score" = sco,
                                   "label" = rownames(p[[sco]]),
                                   weights = wts)
        }
    }
    if (tidy && length(display) > 1) {
        p <- do.call("rbind", p)
        rownames(p) <- p$label
    }
    if (length(p) == 1)
        p[[1]]
    else
        p
}

`plot.wcmdscale` <-
    function(x, choices = c(1,2), display = "sites", type = "t", ...)
{
    ordiplot(x, display = display, choices = choices, type = type, ...)
}

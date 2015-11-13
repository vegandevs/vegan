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
        cat("(Showed only", PRINLIM, "of all", neig, "eigenvalues)\n")
    wvar <- var(x$weights)
    wlen <- length(x$weights)
    cat("\nWeights:")
    if (wvar < 1e-6)
        cat(" Constant\n")
    else {
        cat("\n")
        print(zapsmall(x$weights[1 : min(wlen, PRINLIM)], digits = digits, ...))
        if (wlen > PRINLIM)
            cat("(Showed only", PRINLIM, "of all", wlen, "weights)\n")
    }
    cat("\n")
    invisible(x)
}

`scores.wcmdscale` <-
    function(x, choices = NA, ...)
{
    if (any(is.na(choices)))
        x$points
    else {
        choices <- choices[choices <= NCOL(x$points)]
        x$points[, choices, drop = FALSE]
    }
}

`plot.wcmdscale` <-
    function(x, choices = c(1,2), type = "t", ...)
{
    ordiplot(x, display = "sites", choices = choices, type = type, ...)
}

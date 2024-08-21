`print.cca` <-
    function (x, digits = max(3, getOption("digits") - 3), ...)
{
    ord_obj <- deparse(substitute(x))
    msg_w <- 0.95 * getOption("width")
    writeLines(strwrap(pasteCall(x$call), width = msg_w))
    cat("\n")
    if (!is.null(x$CA$imaginary.chi))
        totchi <- x$tot.chi - x$CA$imaginary.chi
    else
        totchi <- x$tot.chi
    chi <- c(totchi, x$pCCA$tot.chi, x$CCA$tot.chi, x$CA$tot.chi)
    props <- chi/chi[1]
    rnk <- c(NA, x$pCCA$rank, x$CCA$rank, x$CA$rank)
    ## report no. of real axes in dbrda if any negative eigenvalues
    if (inherits(x, "dbrda") &&
        (!is.null(x$CCA) && x$CCA$poseig < x$CCA$qrank ||
             !is.null(x$CA) && x$CA$poseig < x$CA$rank))
        poseig <- c(NA, if (!is.null(x$pCCA)) NA, x$CCA$poseig, x$CA$poseig)
    else
        poseig <- NULL
    tbl <- cbind(chi, props, rnk, poseig)
    if (!is.null(x$CA$imaginary.chi))
        tbl <- rbind(c(x$tot.chi, NA, NA),
                     tbl,
                     c(x$CA$imaginary.chi, NA, NA))
    colnames(tbl) <- c("Inertia", "Proportion", "Rank",
                       if (!is.null(poseig)) "RealDims")
    rn <- c("Total", "RealTotal", "Conditional", "Constrained", "Unconstrained",
            "Imaginary")
    rownames(tbl) <- rn[c(TRUE,
                          !is.null(x$CA$imaginary.chi),
                          !is.null(x$pCCA),
                          !is.null(x$CCA),
                          !is.null(x$CA),
                          !is.null(x$CA$imaginary.chi))]
    ## Remove "Proportion" if only one component
    if (is.null(x$CCA) && is.null(x$pCCA))
        tbl <- tbl[,-2]
    ## 'cs' columns before "Rank" are non-integer
    cs <- which(colnames(tbl) == "Rank") - 1
    writeLines("-- Model Summary --\n")
    printCoefmat(tbl, digits = digits, na.print = "", cs.ind = seq_len(cs))
    writeLines(strwrap(paste("Inertia is", x$inertia), width = msg_w,
        initial = "\n"))
    sp.na <- if (is.null(x$CCA)) attr(x$CA$v, "na.action")
    else attr(x$CCA$v, "na.action")
    # print any notices
    if (any(!is.null(x$vdata), !is.null(x$CCA$alias), x$CA$rank < 1,
        !is.null(x$na.action), !is.null(sp.na))) {
        writeLines("\n-- Note --")
    }
    ## data used for species scores in db ordination
    if (!is.null(x$vdata)) {
        writeLines(
            strwrap(paste("Species scores projected from", sQuote(x$vdata)),
                width = msg_w, initial = "\n")
        )
    }
    # notify if any terms are linearly dependent (aliased)
    if (!is.null(x$CCA$alias)) {
        vif_msg <- sQuote(paste0("vif.cca(", ord_obj, ")"))
        aliased <- paste(sQuote(alias(x, names.only = TRUE)), collapse = ", ")
        msg <- paste("Some constraints or conditions were aliased because they were redundant.",
        "This can happen if terms are linearly dependent (collinear):", aliased)
        writeLines(strwrap(msg, width = msg_w, initial = "\n"))
    }
    if (x$CA$rank < 1) {
        msg <- "The model is overfitted with no unconstrained (residual) component."
        writeLines(strwrap(msg, width = msg_w, initial = "\n"))
    }
    ## Report removed observations and species
    if (!is.null(x$na.action)) {
        writeLines(strwrap(naprint(x$na.action), width = msg_w, initial = "\n"))
    }
    if (!is.null(sp.na)) {
        writeLines(strwrap(paste(length(sp.na), "species",
            ifelse(length(sp.na) == 1, "(variable)", "(variables)"),
            "deleted due to missingness."), width = msg_w, initial = "\n"))
    }
    writeLines("\n-- Eigenvalues --")
    if (!is.null(x$CCA) && x$CCA$rank > 0) {
        cat("\nEigenvalues for constrained axes:\n")
        print(zapsmall(x$CCA$eig, digits = digits), ...)
    }
    if (!is.null(x$CA) && x$CA$rank > 0) {
        ax.lim <- 8
        ax.trig <- 16
        cat("\nEigenvalues for unconstrained axes:\n")
        if (x$CA$rank > ax.trig) {
            print(zapsmall(x$CA$eig[1:ax.lim], digits = digits), ...)
            cat("(Showing", ax.lim, "of", x$CA$rank,
                "unconstrained eigenvalues)\n")
        }
        else print(zapsmall(x$CA$eig, digits = digits), ...)
    }
    cat("\n")
    if (inherits(x, c("capscale", "dbrda"))) {
        if (!is.null(x$metaMDSdist))
            cat("metaMDSdist transformed data:", x$metaMDSdist, "\n\n")
        if (!is.null(x$ac))
            cat("Constant added to distances:", x$ac, "\n\n")
    }
    invisible(x)
}

### package klaR has another function called rda(), and to avoid using
### klar:::print.rda instead of delegating to vegan:::print.cca we
### define here a (redundant) print.rda

`print.rda` <-
    function(x, ...)
{
    ## not vegan rda?
    if (!("CA" %in% names(x)))
        stop(gettextf("%s is not a vegan rda object",
                      sQuote(deparse(substitute(x)))))
    NextMethod("print", x, ...)
}

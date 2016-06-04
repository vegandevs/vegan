`print.cca` <-
    function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    if (inherits(x, "pcaiv")) {
        warning("this is an ade4 object which vegan cannot handle")
        x <- ade2vegancca(x)
    }
    writeLines(strwrap(pasteCall(x$call)))
    cat("\n")
    chi <- c(x$tot.chi, x$pCCA$tot.chi, x$CCA$tot.chi, x$CA$tot.chi)
    props <- chi/chi[1]
    rnk <- c(NA, x$pCCA$rank, x$CCA$rank, x$CA$rank)
    ## handle negative eigenvalues of capscale
    if (!is.null(x$CA$imaginary.chi)) 
        rchi <- c(x$real.tot.chi, x$pCCA$real.tot.chi,
                  x$CCA$real.tot.chi, x$CA$real.tot.chi)
    else
        rchi <- NULL
    ## report no. of real axes in dbrda if any negative eigenvalues
    if (inherits(x, "dbrda") &&
        (!is.null(x$CCA) && x$CCA$poseig < x$CCA$qrank ||
             !is.null(x$CA) && x$CA$poseig < x$CA$rank))
        poseig <- c(NA, if (!is.null(x$pCCA)) NA, x$CCA$poseig, x$CA$poseig)
    else
        poseig <- NULL
    tbl <- cbind(chi, props, rchi, rnk, poseig)
    if (!is.null(rchi))
        tbl <- rbind(tbl, c(NA, NA, x$CA$imaginary.chi,
                            x$CA$imaginary.rank))
    colnames(tbl) <- c("Inertia", "Proportion",
                       if(!is.null(rchi)) "Eigenvals", "Rank",
                       if (!is.null(poseig)) "RealDims")
    rn <- c("Total", "Conditional", "Constrained", "Unconstrained",
            "Imaginary")
    rownames(tbl) <- rn[c(TRUE,!is.null(x$pCCA),
                          !is.null(x$CCA),  !is.null(x$CA),
                          !is.null(x$CA$imaginary.chi))]
    ## Remove "Proportion" if only one component
    if (is.null(x$CCA) && is.null(x$pCCA))
        tbl <- tbl[,-2]
    ## 'cs' columns before "Rank" are non-integer
    cs <- which(colnames(tbl) == "Rank") - 1
    printCoefmat(tbl, digits = digits, na.print = "", cs.ind = seq_len(cs))
    cat("Inertia is", x$inertia, "\n")
    if (!is.null(x$CCA$alias))
        cat("Some constraints were aliased because they were collinear (redundant)\n")
    ## Report removed observations and species
    if (!is.null(x$na.action))
        cat(naprint(x$na.action), "\n")
    sp.na <- if (is.null(x$CCA)) attr(x$CA$v, "na.action")
    else attr(x$CCA$v, "na.action")
    if (!is.null(sp.na))
        cat(length(sp.na), "species",
            ifelse(length(sp.na)==1, "(variable)", "(variables)"),
            "deleted due to missingness\n")
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
            cat("(Showed only", ax.lim, "of all", x$CA$rank, 
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

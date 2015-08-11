`print.cca` <-
    function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    if (inherits(x, "pcaiv")) {
        warning("this is an ade4 object which vegan cannot handle")
        x <- ade2vegancca(x)
    }
    writeLines(strwrap(pasteCall(x$call)))
    cat("\n")
    chi <- c(x$tot.chi, if (!is.null(x$CA$imaginary.chi)) x$tot.chi - x$CA$imaginary.chi,
                            x$pCCA$tot.chi, x$CCA$tot.chi, x$CA$tot.chi,
             x$CA$imaginary.chi)
    ## Proportions of inertia only for Real dimensions in capscale
    if (is.null(x$CA$imaginary.chi))
        props <- chi/chi[1]
    else
        props <- c(NA, chi[-c(1, length(chi))]/chi[2], NA)
    rnk <- c(NA, if (!is.null(x$CA$imaginary.rank)) NA, x$pCCA$rank, x$CCA$rank, x$CA$rank,
             x$CA$imaginary.rank)
    tbl <- cbind(chi, props, rnk)
    colnames(tbl) <- c("Inertia", "Proportion", "Rank")
    rn <- c("Total", "Real Total",  "Conditional", "Constrained", "Unconstrained",
            "Imaginary")
    rownames(tbl) <- rn[c(TRUE, !is.null(x$CA$imaginary.chi), !is.null(x$pCCA),
                          !is.null(x$CCA),  !is.null(x$CA),
                          !is.null(x$CA$imaginary.chi))]
    ## Remove "Proportion" if only one component
    if (is.null(x$CCA) && is.null(x$pCCA))
        tbl <- tbl[,-2]
    printCoefmat(tbl, digits = digits, na.print = "")
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
    invisible(x)
}

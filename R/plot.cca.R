`plot.cca` <- function (x, choices = c(1, 2), display = c("sp", "wa", "cn"),
                        scaling = "species", type, xlim, ylim, const,
                        correlation = FALSE, hill = FALSE,
                        optimize = FALSE, arrows = FALSE,
                        spe.par = list(), sit.par = list(), con.par = list(),
                        bip.par = list(), cen.par = list(), reg.par = list(),
                        ...)
{
    TYPES <- c("text", "points", "none")
    ## take care that bp arrows are also returned if only cn given
    if (any(display %in% c("c","cn")))
        display <- c(display, "bp")
    g <- scores(x, choices, display, scaling, const, correlation = correlation,
                hill = hill, tidy = FALSE, droplist = FALSE)
    if (length(g) == 0 || all(is.na(g)))
      stop("nothing to plot: requested scores do not exist")
    if (!is.list(g)) # never! but doesn't harm either
        g <- list(default = g)
    ## Take care that there are names
    for (i in seq_along(g)) {
        if (length(dim(g[[i]])) > 1)
            rownames(g[[i]]) <- rownames(g[[i]], do.NULL = FALSE,
                                         prefix = substr(names(g)[i], 1, 3))
    }
    if (!is.null(g$centroids)) {
        if (is.null(g$biplot)) # should never be null if g$centroids exist
            g$biplot <- scores(x, choices, "bp", scaling)
        bipnam <- rownames(g$biplot)
        cntnam <- rownames(g$centroids)
        g$biplot <- g$biplot[!(bipnam %in% cntnam), , drop = FALSE]
        if (nrow(g$biplot) == 0)
            g$biplot <- NULL
    }
    if (missing(type)) {
        nitlimit <- 80
        nit <- max(nrow(g$spe), nrow(g$sit), nrow(g$con), nrow(g$cen), 0)
        if (nit > nitlimit)
            type <- "points"
        else type <- "text"
    }
    else type <- match.arg(type, TYPES)
    ## use linestack (and exit) if only one axis was chosen, and
    ## display includes row or column scores
    if (length(choices) == 1) {
        ## Only one set of scores: plot them
        if (length(g) == 1)
            pl <- linestack(g[[1]], ...)
        ## The order of scores is species, sites, constraints, biplot,
        ## centroids: plot two first in this order, but species scores
        ## on the left
        else {
            hasSpec <- names(g)[1] == "species"
            ylim <- range(c(g[[1]], g[[2]]), na.rm = TRUE)
            pl <- linestack(g[[1]], ylim = ylim,
                            side = ifelse(hasSpec, "left", "right"), ...)
            linestack(g[[2]], ylim = ylim,
                      side = ifelse(hasSpec, "right","left"),
                      add = TRUE, ...)
            }
        return(invisible(pl))
    }
    if (missing(xlim)) {
        xlim <- range(0, g$species[, 1], g$sites[, 1], g$constraints[, 1],
                      g$biplot[, 1], g$regression[,1],
                      if (length(g$centroids) > 0 && all(is.na(g$centroids))) NA else g$centroids[, 1],
                      na.rm = TRUE)
    }
    if (!any(is.finite(xlim)))
        stop("no finite scores to plot")
    if (missing(ylim)) {
        ylim <- range(0, g$species[, 2], g$sites[, 2], g$constraints[, 2],
                      g$biplot[, 2], g$regression[,2],
                      if (length(g$centroids) > 0 && all(is.na(g$centroids))) NA else g$centroids[, 2],
                      na.rm = TRUE)
    }
    ordiArgAbsorber(g[[1]], xlim = xlim, ylim = ylim, type = "n", asp = 1,
         FUN = plot, ...)
    abline(h = 0, lty = 3)
    abline(v = 0, lty = 3)
    ## set up lists for graphical parameters
    GlobalPar <- list("type" = type, "optimize" = optimize, "arrows" = arrows)
    dots <- match.call(expand.dots = FALSE)$...
    if (!is.null(dots)) {
        GlobalPar <- modifyList(GlobalPar, dots)
    }
    ## Default graphical parameters
    defParText <- list("species" = list("col" = 2, "cex" = 0.7),
                       "sites" = list("cex" = 0.7),
                       "constraints" = list("col" = "darkgreen", "cex" = 0.7),
                       "biplot" = list("col" = "blue", "cex" = 1.0),
                       "regression" = list("col" = "purple4", "cex" = 1.0),
                       "centroids" = list("col" = "blue", "cex" = 1.0))
    defParPoints <- list("species" = list("col" = 2, "cex" = 0.7, "pch" = "+"),
                         "sites" = list("cex" = 0.7, pch = 1),
                         "constraints" = list("col" = "darkgreen", "cex" = 0.7,
                                              "pch" = 2),
                         "biplot" = list("col" = "blue"),
                         "regression" = list("col" = "purple4"),
                         "centroids" = list("col" = "blue", "pch" = "x",
                                            "cex" = 1.0))
    UserPar <- list("species" = spe.par,
                    "sites" = sit.par,
                    "constraints" = con.par,
                    "biplot" = bip.par,
                    "regression" = reg.par,
                    "centroids" = cen.par)
    ## Plot each score in g. type needs a bit more juggling since it
    ## can be set either globally or for a single score type as user
    ## parameter, and it is not a universal graphical parameter and
    ## must be removed from the final call.
    class(g) <- "ordiplot"
    if (type == "none") return(invisible(g))
    for (kind in names(g)) {
        score <- if (!is.null(UserPar[[kind]]$type))
                     UserPar[[kind]]$type
                 else type
        score <- match.arg(score, TYPES)
        if (score == "none") next
        par <- switch(score,
                      "text" = defParText[[kind]],
                      "points" = defParPoints[[kind]])
        par <- modifyList(par, GlobalPar)
        if (!is.null(UserPar[[kind]]))
            par <- modifyList(par, UserPar[[kind]])
        ## sanitize par combinations
        if (score == "points") # points cannot be optimized
            par <- modifyList(par, list(optimize = NULL))
        else if (score == "text") {
            if (isTRUE(par$optimize)) {
                if (isTRUE(par$arrows) || kind %in% c("biplot", "regression"))
                    message("'optimize = TRUE' and arrows do not mix nicely")
                if (is.null(par$pch)) # optimize=TRUE needs points
                    par <- modifyList(par, list(pch = defParPoints[[kind]]$pch))
            }
        }
        ## add arguments for text/points.ordiplot, remove type
        par <- modifyList(par, list("x" = g, "what" = kind, "type" = NULL))
        do.call(score, par)
    }
    invisible(g)
}

## vegan::plot.rda needed because klaR::plot.rda would be used
## instead if klaR package is loaded

`plot.rda`<-
    function(x, ...)
{
    ## not vegan rda?
    if (!("CA" %in% names(x)))
        stop(gettextf("%s is not a vegan rda object",
                      sQuote(deparse(substitute(x)))))
    NextMethod() # no need to pass ..., happens automagically
}

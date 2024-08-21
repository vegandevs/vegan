`ordiplot` <-
    function (ord, choices = c(1, 2), type = "points", display, xlim,
              ylim, cex = 0.7, ...)
{
    ## local functions to absorb non-par arguments of plot.default
    localPoints <- function(..., log, frame.plot, panel.first,
                            panel.last, axes) points(...)
    localText <- function(..., log, frame.plot, panel.first,
                          panel.last, axes) text(...)
    if (inherits(ord, "decorana") || inherits(ord, "cca")) {
        if (missing(display))
            out <- plot(ord, choices = choices, type = type, xlim = xlim,
                        ylim = ylim, cex = cex, ...)
        else out <- plot(ord, choices = choices, type = type, display = display,
                         xlim = xlim, ylim = ylim, cex = cex, ...)
    }
    else {
        type <- match.arg(type, c("points", "text", "none"))
        ## Matching displays could be done better (see
        ## ordipointlabel), but this may not be yet broken, so...
        dplays <- c("sites", "species")
        if (missing(display))
            display <- dplays
        else
            display <- match.arg(display, dplays, several.ok = TRUE)
        X <- Y <- NULL
        if ("sites" %in% display)
            X <- scores(ord, choices = choices, display = "sites")
        if ("species" %in% display) {
            op <- options(show.error.messages = FALSE)
            Y <- try(scores(ord, choices = choices, display = "species"))
            options(op)
            if (!length(Y) || inherits(Y, "try-error")) {
                message("species scores not available")
                Y <- NULL
            }
            else if (!is.null(X) && NROW(X) == NROW(Y) &&
                     isTRUE(all.equal.numeric(X, Y,
                                              check.attributes = FALSE))) {
                Y <- NULL
                message("species scores not available")
            }
        }
        if (is.null(X) && is.null(Y) || !length(X) && !length(Y))
            stop("no scores found: nothing to plot")
        ## Use linestack and exit if there is only one dimension
        if (NCOL(X) == 1 && NCOL(Y) == 1) {
            pl <- linestack(X, ylim = range(c(X,Y), na.rm=TRUE), cex = cex, ...)
            if (!is.null(Y))
                linestack(Y, side = "left", add = TRUE, cex = cex, ...)
            return(invisible(pl))
        }
        tmp <- apply(rbind(X, Y), 2, range, na.rm=TRUE)
        if (missing(xlim))
            xlim <- tmp[, 1]
        if (missing(ylim))
            ylim <- tmp[, 2]
        plot(tmp, xlim = xlim, ylim = ylim, asp = 1, type = "n",
             ...)
        if (type == "points") {
            if (!is.null(X))
                localPoints(X, pch = 1, col = 1, cex = cex, ...)
            if (!is.null(Y))
                localPoints(Y, pch = "+", col = "red", cex = cex, ...)
        }
        if (type == "text") {
            if (!is.null(X)) {
                labs <- rownames(X)
                if (is.null(labs)) {
                    warning("type='t', but no names available: using x1...")
                    labs <- paste0("x", as.character(seq_len(nrow(X))))
                }
                localText(X, labels = labs, col = 1, cex = cex, ...)
            }
            if (!is.null(Y)) {
                labs <- rownames(Y)
                if (is.null(labs)) {
                    warning("type='t', but no names available: using y1...")
                    labs <- paste0("y", as.character(seq_len(nrow(Y))))
                }
                localText(Y, labels = labs, col = "red", cex = cex, ...)
            }
        }
        out <- list(sites = X, species = Y)
    }
    class(out) <- "ordiplot"
    invisible(out)
}

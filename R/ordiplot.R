`ordiplot` <-
    function (ord, choices = c(1, 2), type = "points", display, xlim, 
              ylim, ...) 
{
    if (!is.null(attr(ord, "class")) && (class(ord) == "decorana" || 
                                         any(class(ord) == "cca"))) {
        if (missing(display)) 
            out <- plot(ord, choices, type = type, xlim = xlim, 
                        ylim = ylim, ...)
        else out <- plot(ord, choices, type = type, display = display, 
                         xlim = xlim, ylim = ylim, ...)
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
            options(show.error.messages = FALSE)
            Y <- try(scores(ord, choices = choices, display = "species"))
            options(show.error.messages = TRUE)
            if (inherits(Y, "try-error")) {
                warning("Species scores not available")
                Y <- NULL
            }
            else if (!is.null(X) && NROW(X) == NROW(Y) &&
                     isTRUE(all.equal.numeric(X, Y,
                                              check.attributes = FALSE))) {
                Y <- NULL
                warning("Species scores not available")
            }
        }
        ## Use linestack and exit if there is only one variable
        if (NCOL(X) == 1 && NCOL(Y) == 1) {
            pl <- linestack(X, ylim = range(c(X,Y), na.rm=TRUE), ...)
            if (!is.null(Y))
                linestack(Y, side = "left", add = TRUE, ...)
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
                points(X, pch = 1, col = 1, cex = 0.7, ...)
            if (!is.null(Y)) 
                points(Y, pch = "+", col = "red", cex = 0.7, 
                       ...)
        }
        if (type == "text") {
            if (!is.null(X)) 
                text(X, labels = rownames(X), col = 1, cex = 0.7, 
                     ...)
            if (!is.null(Y)) 
                text(Y, labels = rownames(Y), col = "red", cex = 0.7, 
                     ...)
        }
        out <- list(sites = X, species = Y)
    }
    class(out) <- "ordiplot"
    invisible(out)
}

"plot.decorana" <-
    function (x, choices = c(1, 2), origin = TRUE, display = c("both", 
                                                   "sites", "species", "none"), cex = 0.8, cols = c(1, 2), type, xlim, ylim, 
              ...) 
{
    display <- match.arg(display)
    sites <- x$rproj
    specs <- x$cproj
    if (missing(type)) {
        nitlimit <- 120
        nit <- 0
        if (display == "sites" || display == "both") 
            nit <- nit + nrow(sites)
        if (display == "species" || display == "both") 
            nit <- nit + nrow(specs)
        if (nit > nitlimit) 
            type <- "points"
        else type <- "text"
    }
    else type <- match.arg(type, c("text", "points", "none"))
    if (origin) {
        sites <- sweep(x$rproj, 2, x$origin, "-")
        specs <- sweep(x$cproj, 2, x$origin, "-")
    }
    sitnam <- rownames(x$rproj)
    spenam <- rownames(x$cproj)
    sites <- sites[, choices]
    specs <- specs[, choices]
    ## Use linestack if only one dim was specified (and exit)
    if (NCOL(sites) == 1 && NCOL(specs) == 1) {
        pl <- linestack(sites,
                        ylim = range(c(sites, specs), na.rm = TRUE), ...)
        linestack(specs, side = "left", add = TRUE, ...)
        return(invisible(pl))
    }
    sp.x <- range(specs[, 1])
    sp.y <- range(specs[, 2])
    st.x <- range(sites[, 1])
    st.y <- range(sites[, 2])
    switch(display, both = {
        if (missing(xlim)) xlim <- range(sp.x, st.x)
        if (missing(ylim)) ylim <- range(sp.y, st.y)
    }, sites = {
        if (missing(xlim)) xlim <- st.x
        if (missing(ylim)) ylim <- st.y
    }, species = {
        if (missing(xlim)) xlim <- sp.x
        if (missing(ylim)) ylim <- sp.y
    }, none = {
        if (missing(xlim)) xlim <- range(sp.x, st.x)
        if (missing(ylim)) ylim <- range(sp.y, st.y)
    })
    plot(sites, type = "n", xlim = xlim, ylim = ylim, asp = 1, 
         ...)
    if (origin) {
        abline(h = 0, lty = 3)
        abline(v = 0, lty = 3)
    }
    else {
        abline(h = x$origin[choices[2]], lty = 3)
        abline(v = x$origin[choices[1]], lty = 3)
    }
    if (type != "none" && (display == "both" || display == "sites")) {
        if (type == "text" && !is.null(sitnam)) 
            text(sites, sitnam, cex = cex, col = cols[1])
        else points(sites, cex = cex, col = cols[1])
    }
    if (type != "none" && (display == "both" || display == "species")) {
        if (type == "text" && !is.null(spenam)) 
            text(specs, spenam, cex = cex, col = cols[2])
        else points(specs, pch = "+", cex = cex, col = cols[2])
    }
    out <- list(sites = sites, species = specs)
    class(out) <- "ordiplot"
    invisible(out)
}

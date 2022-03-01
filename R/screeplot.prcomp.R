`screeplot.prcomp` <-
    function(x, bstick = FALSE, type = c("barplot", "lines"),
             npcs = min(10, length(x$sdev)), ptype = "o", bst.col = "red",
             bst.lty = "solid", xlab = "Component", ylab = "Inertia",
             main = deparse(substitute(x)), legend = bstick, ...)
{
    type <- match.arg(type)
    eig.vals <- x$sdev^2
    ## fix-up names on eig.vals
    names(eig.vals) <- dimnames(x$rotation)[[2]]
    ncomps <- length(eig.vals)
    if(npcs > ncomps)
        npcs <- ncomps
    comps <- seq(len=npcs)
    if(bstick) {
        ord.bstick <- bstick(x)
        ylims <- range(eig.vals[comps], ord.bstick[comps])
    } else {
        ylims <- range(eig.vals)
    }
    if(type=="barplot") {
        ## barplot looks weird if 0 not included
        ylims <- range(0, ylims)
        mids <- barplot(eig.vals[comps], names = names(eig.vals[comps]),
                        main = main, ylab = ylab, ylim = ylims, ...)
    } else {
        plot(comps, eig.vals[comps], type = ptype, axes = FALSE, main = main,
             xlab = xlab, ylab = ylab, ...)
        axis(2)
        axis(1, at = comps, labels = names(eig.vals[comps]))
        mids <-  comps
        box()
    }
    if(bstick) {
        dot.args <- list(...)
        dot.nams <- names(dot.args)
        pch <- if("pch" %in% dot.nams)
            dot.args$pch
        else
            par("pch")
        lines(mids, ord.bstick[comps], type = ptype, col = bst.col,
              lty = bst.lty, pch = pch)
        if(legend) {
            col <- if("col" %in% dot.nams)
                dot.args$col
            else
                par("col")
            lty <- if("lty" %in% dot.nams)
                dot.args$lty
            else
                par("lty")
            if(type == "lines") {
                legend("topright",
                       legend = c("Ordination","Broken Stick"),
                       bty = "n", col = c(col, bst.col),
                       lty = c(lty, bst.lty),
                       pch = pch)
            } else {
                legend("topright",
                       legend = "Broken Stick", bty = "n",
                       col = bst.col, lty = bst.lty, pch = pch)
            }
        }
    }
    invisible(xy.coords(x = mids, y = eig.vals[comps]))
}

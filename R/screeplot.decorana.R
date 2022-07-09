`screeplot.decorana` <-
    function(x, bstick = FALSE, type = c("barplot", "lines"),  npcs = 4,
             ptype = "o", bst.col = "red", bst.lty = "solid",
             xlab = "Component", ylab = "Inertia",
             main = deparse(substitute(x)), legend = bstick, ...)
{
    eig.vals <- if (x$ira == 1) x$evals else x$evals.ortho
    comps <- seq(len=npcs)
    type <- match.arg(type)
    if (bstick) {
        ord.bstick <- bstick(x)
        ylims <- range(0, eig.vals[comps], ord.bstick[comps])
    } else {
        ylims <- range(c(0, eig.vals))
    }
    if(type=="barplot") {
        mids <- barplot(eig.vals[comps], names = names(eig.vals[comps]),
                        main = main, ylab = ylab, ylim = ylims, ...)
    } else {
        plot(eig.vals[comps], type = ptype, axes = FALSE,
             xlab = xlab, ylab = ylab, main = main, ylim = ylims, ...)
        axis(2)
        axis(1, at = comps, labels = names(eig.vals[comps]))
        mids <- comps
        box()
    }
    if (bstick) {
        dot.args <- list(...)
        dot.nams <- names(dot.args)
        pch <- if ("pch" %in% dot.nams)
                   dot.args$pch
               else
                   par("pch")
        lines(mids, ord.bstick[comps], type = ptype, col = bst.col,
              lty = bst.lty, pch = pch)
        if (legend) {
            col <- if ("col" %in% dot.nams)
                       dot.args$col
                   else
                       par("col")
            lty <- if ("lty" %in% dot.nams)
                       dot.args$lty
                   else
                       par("lty")
            if (type == "lines") {
                legend("topright", legend = c("Ordination", "Broken Stick"),
                       bty = "n", col = c(col, bst.col),
                       lty = c(lty, bst.lty), pch = pch)
            }
            else {
                legend("topright", legend = "Broken Stick", bty = "n",
                       col = bst.col, lty = bst.lty, pch = pch)
            }
        }
    }
    invisible(xy.coords(x = mids, y = eig.vals[comps]))
    }

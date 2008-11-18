`screeplot.decorana` <-
    function(x, bstick = FALSE, type = c("barplot", "lines"),  npcs = 4,
             ptype = "o", bst.col = "red", bst.lty = "solid",
             xlab = "Component", ylab = "Inertia",
             main = deparse(substitute(x)), ...)
{
    if (bstick)
        warning("'bstick' not available for 'decorana'")
    eig.vals <- x$evals
    ncomps <- length(eig.vals)
    comps <- seq(len=npcs)
    type <- match.arg(type)
    if(type=="barplot") {
        mids <- barplot(eig.vals[comps], names = names(eig.vals[comps]),
                        main = main, ylab = ylab, ...)
    } else {
        plot(eig.vals[comps], type = ptype, axes = FALSE,
             xlab = xlab, ylab = ylab, main = main, ...)
        axis(2)
        axis(1, at = comps, labels = names(eig.vals[comps]))
        mids <- comps
        box()
    }
    invisible(xy.coords(x = mids, y = eig.vals[comps]))
}

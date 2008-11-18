`screeplot.prcomp` <-
    function(x, bstick = FALSE, type = c("barplot", "lines"),
             npcs = min(10, length(x$sdev)), ptype = "o", bst.col = "red",
             bst.lty = "solid", xlab = "Component", ylab = "Inertia",
             main = deparse(substitute(x)), ...)
{
    main
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
        lines(mids, ord.bstick[comps], type = ptype, col = bst.col,
              lty = bst.lty)
    }
    invisible(xy.coords(x = mids, y = eig.vals[comps]))
}

`lines.spantree` <-
    function (x, ord, display = "sites", col = 1, ...)
{
    ord <- scores(ord, display = display, ...)
    tree <- x$kid
    ## recycle colours and use a mixture of joined points for line segments
    col <- rep(col, length = nrow(ord))
    col <- col2rgb(col)/255
    ## average colour for pairs of points
    col <- rgb(t(col[,-1] + col[,tree])/2)
    if (x$n > 1)
        ordiArgAbsorber(ord[-1, 1], ord[-1, 2], ord[tree, 1], ord[tree, 2],
                        col = col,
                        FUN = segments, ...)
    invisible()
}

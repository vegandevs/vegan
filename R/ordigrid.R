"ordigrid" <-
    function (ord, levels, replicates, display = "sites", lty=c(1,1), col=c(1,1),
              lwd = c(1,1), ...) 
{
    pts <- scores(ord, display = display, ...)
    npoints <- nrow(pts)
    gr <- gl(levels, replicates, npoints)
    ordisegments(pts, groups = gr,  lty = lty[1], col = col[1],
                 lwd = lwd[1], ...)
    gr <- gl(replicates, 1, npoints)
    ordisegments(pts, groups = gr, lty = lty[2], col = col[2],
                 lwd = lwd[2], ...)
    invisible()
}

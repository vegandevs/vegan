"ordigrid" <-
    function (ord, levels, replicates, display = "sites", ...) 
{
    pts <- scores(ord, display = display, ...)
    npoints <- nrow(pts)
    gr <- gl(levels, replicates, npoints)
    ordisegments(pts, groups = gr, ...)
    gr <- gl(replicates, 1, npoints)
    ordisegments(pts, groups = gr, ...)
    invisible()
}

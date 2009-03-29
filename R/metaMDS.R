`metaMDS` <-
    function (comm, distance = "bray", k = 2, trymax = 20, autotransform = TRUE, 
              noshare = 0.1, wascores = TRUE, expand = TRUE, trace = 1,
              plot = FALSE, previous.best, old.wa = FALSE, ...) 
{
    commname <- deparse(substitute(comm))
    if (inherits(comm, "dist"))
        stop("metaMDS does not accept dist objects, but needs original data frame")
    dis <- metaMDSdist(comm, distance = distance, autotransform = autotransform, 
                       noshare = noshare, trace = trace, commname = commname, 
                       ...)
    if (missing(previous.best)) 
        previous.best <- NULL
    out <- metaMDSiter(dis, k = k, trymax = trymax, trace = trace, 
                       plot = plot, previous.best = previous.best, ...)
    points <- postMDS(out$points, dis, plot = max(0, plot - 1), ...)
    if (is.null(rownames(points))) 
        rownames(points) <- rownames(comm)
    if (wascores) {
        if (!old.wa)
            comm <- eval.parent(parse(text=attr(dis, "commname")))
        wa <- wascores(points, comm, expand = expand)
        attr(wa, "old.wa") <- old.wa
    }
    else
        wa <- NA
    out$points <- points
    out$species <- wa
    out$call <- match.call()
    class(out) <- "metaMDS"
    out
}

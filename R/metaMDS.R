`metaMDS` <-
    function (comm, distance = "bray", k = 2, trymax = 20,
              engine = c("monoMDS", "isoMDS"), 
              autotransform = TRUE, noshare = (engine == "isoMDS"),
              wascores = TRUE, expand = TRUE, trace = 1,
              plot = FALSE, previous.best,  ...) 
{
    engine <- match.arg(engine)
    ## This could be a character vector of length > 1L
    commname <- deparse(substitute(comm), width.cutoff = 500L)
    if (length(commname) > 1L) {
        paste(commname, collapse = "", sep = "")
        ## deparse can add more white space, so cull 2 or more spaces to a single space
        commname <- gsub("[ ]{2,}", " ", commname)
    }
    ## metaMDS was written for community data which should be all
    ## positive. Check this here, and set arguments so that they are
    ## suitable for non-negative data.
    if (any(autotransform, noshare > 0, wascores) && any(comm < 0, na.rm=TRUE)) {
        warning("'comm' has negative data: 'autotransform', 'noshare' and 'wascores' set to FALSE")
        wascores <- FALSE
        autotransform <- FALSE
        noshare <- FALSE
    }
    if (inherits(comm, "dist")) {
        dis <- comm
        if (is.null(attr(dis, "method")))
            attr(dis, "method") <- "user supplied"
        wascores <- FALSE
    } else if (length(dim(comm) == 2) && ncol(comm) == nrow(comm) &&
                all(comm == t(comm))) {
        dis <- as.dist(comm)
        attr(dis, "method") <- "user supplied"
        wascores <- FALSE
    } else {
        dis <- metaMDSdist(comm, distance = distance,
                           autotransform = autotransform, 
                           noshare = noshare, trace = trace,
                           commname = commname, ...)
    }
    if (missing(previous.best)) 
        previous.best <- NULL
    out <- metaMDSiter(dis, k = k, trymax = trymax, trace = trace, 
                       plot = plot, previous.best = previous.best,
                       engine = engine, ...)
    ## Nearly zero stress is usually not a good thing but a symptom of
    ## a problem: you may have insufficient data for NMDS
    if (out$stress < 1e-3) {
        warning("Stress is (nearly) zero - you may have insufficient data")
    }     
    points <- postMDS(out$points, dis, plot = max(0, plot - 1), ...)
    if (is.null(rownames(points))) 
        rownames(points) <- rownames(comm)
    wa <- if (wascores) {
        ## transformed data
        ##comm <- eval.parent(parse(text=attr(dis, "commname")))
        comm <- attr(dis, "comm")
        wascores(points, comm, expand = expand)
    } else {
        NA
    }
    out$points <- points
    out$species <- wa
    out$call <- match.call()
    if (is.null(out$data))
        out$data <- commname
    class(out) <- c("metaMDS", class(out))
    out
}

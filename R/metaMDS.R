`metaMDS` <-
    function (comm, distance = "bray", k = 2, try = 20, trymax = 20,
              engine = monoMDS, autotransform = TRUE, noshare = FALSE,
              wascores = TRUE, expand = TRUE, trace = 1,
              plot = FALSE, previous.best,  ...)
{
    ## take care that try (minimum) is not larger than trymax (maximum)
    if (try > trymax)
        try <- trymax
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
        message("'comm' has negative data: 'autotransform', 'noshare' and 'wascores' set to FALSE")
        wascores <- FALSE
        autotransform <- FALSE
        noshare <- FALSE
    }
    if (inherits(comm, "dist")) {
        dis <- comm
        if (is.null(attr(dis, "method")))
            attr(dis, "method") <- "user supplied"
        wascores <- FALSE
    } else if ((is.matrix(comm) || is.data.frame(comm)) &&
               isSymmetric(unname(as.matrix(comm)))) {
        dis <- as.dist(comm)
        attr(dis, "method") <- "user supplied"
        wascores <- FALSE
    } else {
        if (trace > 2)
            cat(">>> Calculation of dissimilarities\n")
        dis <- metaMDSdist(comm, distance = distance,
                           autotransform = autotransform,
                           noshare = noshare, trace = trace,
                           commname = commname, ...)
    }
    if (missing(previous.best))
        previous.best <- NULL
    if (trace > 2)
        cat(">>> NMDS iterations\n")
    maker <- deparse1(substitute(engine))
    out <- metaMDSiter(dis, k = k, try = try, trymax = trymax,
                       trace = trace,
                       plot = plot, previous.best = previous.best,
                       engine = engine, maker = maker, ...)
    out$engine <- maker
    ## Nearly zero stress is usually not a good thing but a symptom of
    ## a problem: you may have insufficient data for NMDS
    if (out$stress < 1e-3) {
        warning("stress is (nearly) zero: you may have insufficient data")
    }
    if (trace > 2)
        cat(">>> Post-processing NMDS\n")
    points <- postMDS(out$points, dis, plot = max(0, plot - 1), ...)
    ## rescale monoMDS scaling if postMDS scaled 'points'
    if (out$engine == "monoMDS" &&
        !is.null(scl <- attr(points, "internalscaling"))) {
        out$dist <- out$dist/scl
        out$dhat <- out$dhat/scl
    }
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

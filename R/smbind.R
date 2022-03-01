`smbind` <-
    function (object, ..., MARGIN, strict = TRUE)
{
    if (missing(MARGIN))
        stop("MARGIN argument must be specified")
    MARGIN <- as.integer(MARGIN)
    if (length(MARGIN) != 1L)
        stop("MARGIN length must be 1")
    if (!(MARGIN %in% 1L:3L))
        stop("MARGIN value must be in 1:3")

    if (is.list(object)) {
        obj <- object
        if (!missing(...))
            warning("'object' was a list, '...' ignored")
    } else {
        obj <- list(object, ...)
    }
    l <- length(obj)
    if (l < 2L)
        return(obj[[1L]])
    att <- lapply(obj, attributes)
    isSeq <- att[[1L]]$isSeq
    startEq <- endEq <- thinEq <- OKseed <- TRUE
    for (i in 2L:l) {
        ## data must be identical when MARGIN=3
        if (MARGIN == 3L && !identical(att[[1L]][["data"]], att[[i]][["data"]]))
            stop("'data' attributes are not identical")
        ## dimensions need to match except for MARGIN
        if (!identical(att[[1L]][["dim"]][-MARGIN], att[[i]][["dim"]][-MARGIN]))
            stop("dimension mismatch")
        ## method settings need to be set on return object
        ## thus these need to be identical
        for (NAM in c("method", "binary", "isSeq", "mode", "class")) {
            if (!identical(att[[1L]][[NAM]], att[[i]][[NAM]]))
                stop(gettextf("'%s' attributes not identical", NAM))
        }
        ## ts attributes are tricky: evaluate outside of the loop
        for (NAM in c("start", "end", "thin")) {
            if (!identical(att[[1L]][["start"]], att[[i]][["start"]]))
                startEq <- FALSE
            if (!identical(att[[1L]][["end"]], att[[i]][["end"]]))
                endEq <- FALSE
            if (!identical(att[[1L]][["thin"]], att[[i]][["thin"]]))
                thinEq <- FALSE
        }
        ## seed is important when 'data' are the same (MARGIN=3)
        ## but it is up to the user
        ## return value has NULL seed attribute
        if (MARGIN == 3L && identical(att[[1L]][["seed"]], att[[i]][["seed"]])) {
            OKseed <- FALSE
        }
    }
    if (!OKseed)
        warning("identical 'seed' attributes found")
    if (isSeq) {
        outStart <- outEnd <- outThin <- NA
        type <- "none"
        ## if MARGIN != 3
        ##   all match or fail
        ##   when all match: keep ts attributes, type: "strat"
        ## if MARGIN==3
        ##   sequential algorithms need identical ts attributes
        ##   * if parallel (start/end/thin identical): "par"
        ##   --> original start, end, thin, + set chains attr
        ##   * if subsequent (start/end/thin form a sequence): "seq"
        ##   --> calculate start & end, thin same
        ##   * all else: "none"
        ##   --> fail unless strict=FALSE (when start=NA, end=NA, thin=NA)
        if (MARGIN != 3L) {
            if (startEq && endEq && thinEq) {
                type <- "strat"
                outStart <- att[[1L]]$start
                outEnd <- att[[1L]]$end
                outThin <- att[[1L]]$thin
            }
        } else {
            if (startEq && endEq && thinEq) {
                type <- "par"
                outStart <- att[[1L]]$start
                outEnd <- att[[1L]]$end
                outThin <- att[[1L]]$thin
            }
            if (!startEq && !endEq && thinEq) {
                stv <- sapply(att, "[[", "start")
                o <- order(stv)
                att <- att[o]
                obj <- obj[o]
                stv <- sapply(att, "[[", "start")
                env <- sapply(att, "[[", "end")
                thv <- att[[1L]]$thin
                nsv <- sapply(obj, function(z) dim(z)[3L])
                vals <- lapply(1:l, function(i)
                    seq(stv[i], env[i], by=thv))
                OK <- logical(4L)
                if (length(stv) == length(unique(stv)))
                    OK[1L] <- TRUE
                if (length(env) == length(unique(env)))
                    OK[2L] <- TRUE
                if (all(nsv == sapply(vals, length)))
                    OK[3L] <- TRUE
                if (length(seq(stv[1], env[l], by=thv)) == length(unlist(vals)))
                    OK[4L] <- TRUE
                if (all(OK)) {
                    if (all(seq(stv[1], env[l], by=thv) == unlist(vals))) {
                            type <- "seq"
                            outStart <- stv[1]
                            outEnd <- env[l]
                            outThin <- thv
                    }
                }
            }
        }
        if (type == "none") {
            if (strict) {
                stop("incosistent 'start', 'end', 'thin' attributes")
            } else {
                warning("incosistent 'start', 'end', 'thin' attributes")
            }
        }
    }
    ## set final dimensions
    DIM <- att[[1L]]$dim
    DIMs <- sapply(att, function(z) z$dim[MARGIN])
    cDIMs <- cumsum(DIMs)
    DIM[MARGIN] <- cDIMs[l]
    out <- array(NA, dim = DIM)
    ## copy the 1st object
    if (MARGIN == 1L)
        out[1L:dim(obj[[1L]])[1L],,] <- obj[[1L]]
    if (MARGIN == 2L)
        out[,1L:dim(obj[[1L]])[2L],] <- obj[[1L]]
    if (MARGIN == 3L)
        out[,,1L:dim(obj[[1L]])[3L]] <- obj[[1L]]
    ## data attribute will change when MARGIN != 3
    DATA <- att[[1L]]$data
    ## copy 2:l objects and data argument
    for (i in 2L:l) {
        j <- (cDIMs[i - 1L] + 1L):cDIMs[i]
        if (MARGIN == 1L) {
            out[j,,] <- obj[[i]]
            DATA <- rbind(DATA, att[[i]]$data)
        }
        if (MARGIN == 2L) {
            out[,j,] <- obj[[i]]
            DATA <- cbind(DATA, att[[i]]$data)
        }
        if (MARGIN == 3L) {
            out[,,j] <- obj[[i]]
        }
    }
    ## assembling return object
    ratt <- att[[1L]]
    ratt$data <- DATA
    ratt$seed <- NA
    ratt$dim <- DIM
    if (!isSeq)
        ratt$end <- cDIMs[l]
    if (isSeq) {
        ratt$start <- outStart
        ratt$end <- outEnd
        ratt$thin <- outThin
        if (type == "par")
            ratt$chains <- l
    }
    ratt$dimnames[[MARGIN]] <- make.names(unlist(lapply(att, function(z)
        z$dimnames[[MARGIN]])), unique = TRUE)
    attributes(out) <- ratt
    out
}

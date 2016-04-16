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
    if (l < 2)
        return(obj[[1L]])
    att <- lapply(obj, attributes)
    isSeq <- att[[1L]]$isSeq
    OKstart <- OKend <- OKthin <- TRUE
    for (i in 2L:l) {
        ## data must be identical when MARGIN=3
        if (MARGIN == 3L && !identical(att[[1L]][["data"]], att[[i]][["data"]]))
            stop("'data' attributes not identical")
        ## dimensions need to match except for MARGIN
        if (!identical(att[[1L]][["dim"]][-MARGIN], att[[i]][["dim"]][-MARGIN]))
            stop("dimension mismatch")
        ## method settings need to be set on return object
        ## thus these need to be identical
        for (NAM in c("method", "binary", "isSeq", "mode", "class")) {
            if (!identical(att[[1L]][[NAM]], att[[i]][[NAM]]))
                stop("'", NAM, "' attributes not identical")
        }
        ## sequential algorithms need identical ts attributes
        if (isSeq) {
            for (NAM in c("start", "end", "thin")) {
                if (strict && !identical(att[[1L]][[NAM]], att[[i]][[NAM]]))
                    stop("'", NAM, "' attributes not identical")
                if (!strict) {
                    if (!identical(att[[1L]][[NAM]], att[[i]][[NAM]])) {
                        warning("'", NAM, "' attributes not identical")
                        if (NAM == "start")
                            OKstart <- FALSE
                        if (NAM == "end")
                            OKend <- FALSE
                        if (NAM == "thin")
                            OKthin <- FALSE
                    }
                }
            }
        }
        ## seed is important when 'data' are the same (MARGIN=3)
        ## but it is up to the user
        ## return value has NULL seed attribute
        if (MARGIN == 3L &&
            identical(att[[1L]][["seed"]], att[[i]][["seed"]])) {
            warning("identical 'seed' attributes found")
        }
    }
    ## set final dimensions
    DIM <- att[[1L]]$dim
    DIMs <- sapply(att, function(z) z$dim[MARGIN])
    cDIMs <- cumsum(DIMs)
    DIM[MARGIN] <- cDIMs[l]
    out <- array(NA, dim = DIM)
    ## copy in the 1st object
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
    ratt$seed <- NULL
    ratt$dim <- DIM
    if (!isSeq)
        ratt$end <- cDIMs[l]
    if (isSeq) {
        if (!OKstart)
            ratt$start <- NA
        if (!OKend)
            ratt$end <- NA
        if (!OKthin)
            ratt$thin <- NA
    }
    ratt$dimnames[[MARGIN]] <- make.names(unlist(lapply(att, function(z)
        z$dimnames[[MARGIN]])), unique = TRUE)
    attributes(out) <- ratt
    out
}

## permatfull function
`permatfull` <-
function(m, fixedmar="both", shuffle="both",
strata=NULL, mtype="count", times=99, ...)
{
    mtype <- match.arg(mtype, c("prab", "count"))
    shuffle <- match.arg(shuffle, c("ind", "samp", "both"))
    fixedmar <- match.arg(fixedmar, c("none", "rows", "columns", "both"))
    m <- as.matrix(m, force.rownames = TRUE)
    str <- if (is.null(strata))
        1 else as.integer(as.factor(strata)[drop = TRUE])
    levstr <- unique(str)
    nstr <- length(unique(str))
    if (!is.null(strata) && any(table(str) < 2))
        stop("strata should contain at least two observations")
    ALGO <- switch(fixedmar,
        "none" = "r00",
        "rows" = "r0",
        "columns" = "c0",
        "both" = ifelse(mtype=="prab", "quasiswap", "r2dtable"))
    if (mtype=="count") {
        if (fixedmar!="both")
            ALGO <- paste(ALGO, shuffle, sep="_")
    }
    if (is.null(strata)) {
        tmp <- simulate(nullmodel(m, ALGO), nsim=times, ...)
        perm <- vector("list", times)
        for (i in seq_len(times))
            perm[[i]] <- tmp[,,i]
    } else {
        perm <- vector("list", times)
        tmp <- vector("list", length(unique(strata)))
        for (j in seq_len(nstr)) {
            tmp[[j]] <- simulate(nullmodel(m[strata==levstr[j],], ALGO),
                nsim=times, ...)
        }
        for (i in seq_len(times)) {
            perm[[i]] <- array(0, dim(m))
            for (j in seq_len(nstr)) {
                perm[[i]][strata==levstr[j],] <- tmp[[j]][,,i]
            }
        }
    }
    if (fixedmar == "both")
        shuffle <- NA
    if (mtype == "prab")
        m <- ifelse(m > 0, 1, 0)
    out <- list(call=match.call(), orig=m, perm=perm)
    attr(out, "mtype") <- mtype
    attr(out, "ptype") <- "full"
    attr(out, "method") <- ALGO
    attr(out, "fixedmar") <- fixedmar
    attr(out, "times") <- times
    attr(out, "shuffle") <- shuffle
    attr(out, "is.strat") <- !is.null(strata)
    attr(out, "strata") <- str
    attr(out, "burnin") <- NA
    attr(out, "thin") <- NA
    class(out) <- c("permatfull", "permat")
    out
}

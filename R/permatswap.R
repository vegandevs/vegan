## permatswap function
`permatswap` <-
function(m, method="quasiswap", fixedmar="both", shuffle="both", strata=NULL,
mtype="count", times=99, burnin = 0, thin = 1, ...)
{
    mtype <- match.arg(mtype, c("prab", "count"))
    fixedmar <- match.arg(fixedmar, c("none", "rows", "columns", "both"))
    shuffle <- match.arg(shuffle, c("samp", "both"))
    count <- mtype == "count"
    m <- as.matrix(m)
    str <- if (is.null(strata))
        1 else as.integer(as.factor(strata)[drop = TRUE])
    levstr <- unique(str)
    nstr <- length(unique(str))
    if (!is.null(strata) && any(table(str) < 2))
        stop("strata should contain at least two observations")
    ## evaluating algo type
    if (count) {
        method <- match.arg(method, c("swap", "quasiswap", "swsh", "abuswap"))
        if (method == "swap") {
            warning("quantitative swap method may not yield random null models, use only to study its properties")
            isSeq <- TRUE
            if (fixedmar != "both")
                stop("if 'method=\"swap\"', 'fixedmar' must be \"both\"")
        } else {
            if (method == "abuswap") {
                if (fixedmar == "both")
                    stop("if 'method=\"abuswap\"', 'fixedmar' must not be \"both\"")
                direct <- if (fixedmar == "columns")
                    0 else 1
                isSeq <- TRUE
            } else {
                isSeq <- FALSE
                if (method != "swsh" && fixedmar != "both")
                    stop("'fixedmar' must be \"both\"")
            }
        }
        if (method %in% c("swap", "quasiswap"))
            ALGO <- paste(method, "count", sep="_")
        if (method == "abuswap")
            ALGO <- paste(method, substr(fixedmar, 1, 1), sep="_")
        if (method == "swsh") {
            if (fixedmar=="both")
                stop("if 'method=\"swsh\"', 'fixedmar' must not be \"both\"")
            ALGO <- if (fixedmar=="none") {
                paste(method, shuffle, sep="_")
            } else {
                paste(method, shuffle, substr(fixedmar, 1, 1), sep="_")
            }
        }
    } else {
        if (fixedmar != "both")
            stop("if 'mtype=\"prab\"', 'fixedmar' must be \"both\"")
        method <- match.arg(method, c("swap", "quasiswap", "tswap", "backtrack"))
        isSeq <- method != "quasiswap"
        ALGO <- method
    }
    if (is.null(strata)) {
        tmp <- simulate(nullmodel(m, ALGO),
            nsim=times, burnin=burnin, thin=thin, ...)
        perm <- vector("list", times)
        for (i in seq_len(times))
            perm[[i]] <- tmp[,,i]
    } else {
        perm <- vector("list", times)
        tmp <- vector("list", length(unique(strata)))
        for (j in seq_len(nstr)) {
            tmp[[j]] <- simulate(nullmodel(m[strata==levstr[j],], ALGO),
                nsim=times, burnin=burnin, thin=thin, ...)
        }
        for (i in seq_len(times)) {
            perm[[i]] <- array(0, dim(m))
            for (j in seq_len(nstr)) {
                perm[[i]][strata==levstr[j],] <- tmp[[j]][,,i]
            }
        }
    }
    if (mtype == "prab")
        m <- ifelse(m > 0, 1, 0)
    out <- list(call=match.call(), orig=m, perm=perm)
    attr(out, "mtype") <- mtype
    attr(out, "ptype") <- "swap"
    attr(out, "method") <- ALGO
    attr(out, "fixedmar") <- if (method == "swsh") "none" else fixedmar
    attr(out, "times") <- times
    attr(out, "shuffle") <- if (method == "swsh") shuffle else NA
    attr(out, "is.strat") <- !is.null(strata)
    attr(out, "strata") <- str
    attr(out, "burnin") <- burnin
    attr(out, "thin") <- thin
    class(out) <- c("permatswap", "permat")
    out
}

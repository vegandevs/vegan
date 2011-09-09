## permatswap function
`permatswap` <-
function(m, method="quasiswap", fixedmar="both", shuffle="both", strata=NULL,
         mtype="count", times=99, burnin = 0, thin = 1)
{
    mtype <- match.arg(mtype, c("prab", "count"))
    fixedmar <- match.arg(fixedmar, c("rows", "columns", "both"))
    shuffle <- match.arg(shuffle, c("samp", "both"))
    count <- mtype == "count"
    ## evaluating algo type
    if (count) {
        method <- match.arg(method, c("swap", "quasiswap", "swsh", "abuswap"))
        ## warning if swapcount is to be used
        if (method == "swap") {
            warning("quantitative swap method may not yield random null models, use only to study its properties")
            isSeq <- TRUE
        } else {
            if (method == "abuswap") {
                isSeq <- TRUE
            } else {
                isSeq <- FALSE
            }
        }
    } else {
        method <- match.arg(method, c("swap", "quasiswap", "tswap", "backtracking"))
        isSeq <- method != "quasiswap"
    }

    ## sequential algos: might have burnin
    if (isSeq) {
        perm <- vector("list", times)
        ## with burnin
        perm[[1]] <- permatswap1(m, method, fixedmar, shuffle, strata, mtype, burnin+thin)
        for (i in 2:times)
            perm[[i]] <- permatswap1(perm[[i-1]], method, fixedmar, shuffle, strata, mtype, thin)
    ## non-sequential algos: no burnin required
    } else {
        if (burnin > 0)
            warnings("Non sequential algorithm used: 'burnin' argument ignored")
        burnin <- 0
        if (thin != 1)
            warnings("Non sequential algorithm used: 'thin' value set to 1")
        thin <- 1
        perm <- replicate(times, 
            permatswap1(m, method, fixedmar, shuffle, strata, mtype, thin), 
            simplify=FALSE)
    }
    if (mtype == "prab")
        m <- ifelse(m > 0, 1, 0)
    out <- list(call=match.call(), orig=m, perm=perm)
    attr(out, "mtype") <- mtype
    attr(out, "ptype") <- "swap"
    attr(out, "method") <- method
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

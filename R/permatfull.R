## permatfull function
`permatfull` <-
function(m, fixedmar="both", shuffle="both", strata=NULL, mtype="count", times=99)
{
    mtype <- match.arg(mtype, c("prab", "count"))
    shuffle <- match.arg(shuffle, c("ind", "samp", "both"))
    fixedmar <- match.arg(fixedmar, c("none", "rows", "columns", "both"))
    m <- as.matrix(m)

    if (is.null(strata))
        str <- as.factor(rep(1, nrow(m)))
        else str <- as.factor(strata)[drop = TRUE]

    levels(str) <- 1:length(unique(str))
    str <- as.numeric(str)
    nstr <- length(unique(str))
    if (any(tapply(str,list(str),length) == 1))
        stop("strata should contain at least 2 observations")

    perm <- replicate(times, permatfull1(m, fixedmar, shuffle, strata, mtype), simplify=FALSE)
    if (fixedmar == "both")
        shuffle <- NA
    if (mtype == "prab")
        m <- ifelse(m > 0, 1, 0)
    out <- list(call=match.call(), orig=m, perm=perm)
    attr(out, "mtype") <- mtype
    attr(out, "ptype") <- "full"
    attr(out, "method") <- NA
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

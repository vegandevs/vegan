## permatfull function
`permatfull` <-
function(m, fixedmar="both", shuffle="both", strata=NULL, mtype="count", times=99)
{
## internal function
indshuffle <- function(x)
{
   N <- length(x)
   n <- sum(x)
   out <- numeric(N)
   names(out) <- 1:N
   y <- table(sample(1:N, n, replace = TRUE))
   out[names(out) %in% names(y)] <- y
   names(out) <- NULL
   return(out)
}
bothshuffle <- function(x, y=1)
{
    x[x!=0] <- indshuffle(x[x!=0] - y) + y
    return(sample(x))
}
    if (!identical(all.equal(m, round(m)), TRUE))
       stop("function accepts only integers (counts)")
    mtype <- match.arg(mtype, c("prab", "count"))
    shuffle <- match.arg(shuffle, c("ind", "samp", "both"))
    count <- mtype == "count"
    fixedmar <- match.arg(fixedmar, c("none", "rows", "columns", "both"))
    sample.fun <- switch(shuffle,
        "ind"=indshuffle,
        "samp"=sample,
        "both"=bothshuffle)
    m <- as.matrix(m)
    n.row <- nrow(m)
    n.col <- ncol(m)
    if (mtype == "prab") m <- ifelse(m > 0, 1, 0)

    if (is.null(strata))
        str <- as.factor(rep(1, n.row))
        else str <- as.factor(strata)[drop = TRUE]

    levels(str) <- 1:length(unique(str))
    str <- as.numeric(str)
    nstr <- length(unique(str))
    if (any(tapply(str,list(str),length) == 1))
        stop("strata should contain at least 2 observations")
    perm <- list()
    perm[[1]] <- matrix(0, n.row, n.col)
    for (k in 1:times)
        perm[[k]] <- perm[[1]]
    for (j in 1:nstr) {
    id <- which(str == j)
        if (fixedmar == "none")
            for (i in 1:times)
                if (count) perm[[i]][id,] <- matrix(sample.fun(array(m[id,])), length(id), n.col)
                else perm[[i]][id,] <- commsimulator(m[id,], method="r00")
        if (fixedmar == "rows")
            for (i in 1:times)
                if (count) perm[[i]][id,] <- t(apply(m[id,], 1, sample.fun))
                else perm[[i]][id,] <- commsimulator(m[id,], method="r0")
        if (fixedmar == "columns")
            for (i in 1:times)
                if (count) perm[[i]][id,] <- apply(m[id,], 2, sample.fun)
                else perm[[i]][id,] <- commsimulator(m[id,], method="c0")
        if (fixedmar == "both")
            for (i in 1:times)
                if (count) perm[[i]][id,] <- r2dtable(1, apply(m[id,], 1, sum), apply(m[id,], 2, sum))[[1]]
                else perm[[i]][id,] <- commsimulator(m[id,], method="quasiswap")
        }
    if (fixedmar == "both")
        shuffle <- NA
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
    return(out)
}

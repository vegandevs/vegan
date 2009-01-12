## permatswap function
`permatswap` <-
function(m, method="quasiswap", shuffle="both", strata=NULL, mtype="count", times=99, burnin = 10000, thin = 1000)
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
    shuffle <- match.arg(shuffle, c("samp", "both"))
    count <- mtype == "count"
    if (count) {
        method <- match.arg(method, c("swap", "quasiswap", "swsh"))
        ## warning if swapcount is to be used
        if (method == "swap") {
            warning("quantitative swap method may not yield random null models, use only to study its properties")
            isSeq <- TRUE
        } else isSeq <- FALSE
    } else {
        method <- match.arg(method, c("swap", "quasiswap", "tswap", "backtracking"))
        isSeq <- method != "quasiswap"
        }

    m <- as.matrix(m)
    att <- attributes(m)
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
    for (i in 2:times)
        perm[[i]] <- perm[[1]]

    for (j in 1:nstr) {
        id <- which(str == j)
        temp <- m[id,]
        nn.row <- nrow(m[id,])
        nn.col <- ncol(m[id,])
        if (isSeq) {
            if (count) {
                for (k in 1:burnin)
                    temp <- .C("swapcount", m = as.double(temp),
                            as.integer(nn.row), as.integer(nn.col),
                            as.integer(1), PACKAGE = "vegan")$m
            } else
                for (k in 1:burnin)
                    temp <- commsimulator(temp, method=method)
            for (i in 1:times) {
                if (count) {
                    perm[[i]][id,] <- .C("swapcount",
                                    m = as.double(temp),
                                    as.integer(nn.row),
                                    as.integer(nn.col),
                                    as.integer(thin),
                                    PACKAGE = "vegan")$m
	           } else perm[[i]][id,] <- commsimulator(temp, method=method, thin=thin)
            temp <- perm[[i]][id,]
            } # for i end
        } else {
            if (method != "swsh") {
                r2tabs <- r2dtable(times, rowSums(m[id,]), colSums(m[id,]))
                } else tempPos <- temp[temp > 0]
            for (i in 1:times) {
                if (count) {
                    if (method != "swsh") {
                        ms <- sum(m[id,] > 0)
                        tmp <- r2tabs[[i]]
                        ## if fills are equal, no need to restore fill
                        if (sum(tmp > 0) != ms) {
                            tmp <- .C("rswapcount",
                                        m = as.double(tmp),
                                        as.integer(nn.row),
                                        as.integer(nn.col),
                                        as.integer(ms),
                                        PACKAGE="vegan")$m
                        }
                        perm[[i]][id,] <- matrix(tmp, nrow(perm[[i]][id,]), ncol(perm[[i]][id,]))
                    } else { # method == "swsh"
                        tmp <- commsimulator(temp, method="quasiswap")
                        if (shuffle == "samp") {
                            tmp[tmp > 0] <- sample(tempPos)
                        } else tmp[tmp > 0] <- bothshuffle(tempPos)
                        perm[[i]][id,] <- tmp
                    }
                } else perm[[i]][id,] <- commsimulator(temp, method=method)
            }
            thin <- burnin <- 0
        }
    } # for j end
    out <- list(call=match.call(), orig=m, perm=perm)
    attr(out, "mtype") <- mtype
    attr(out, "ptype") <- "swap"
    attr(out, "method") <- method
    attr(out, "fixedmar") <- if (method == "swsh") "none" else "both"
    attr(out, "times") <- times
    attr(out, "shuffle") <- if (method == "swsh") shuffle else NA
    attr(out, "is.strat") <- !is.null(strata)
    attr(out, "strata") <- str
    attr(out, "burnin") <- burnin
    attr(out, "thin") <- thin
    class(out) <- c("permat", "list")
    return(out)
}

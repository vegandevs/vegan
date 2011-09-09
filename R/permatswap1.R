## permatswap1 returns 1 permuted matrix
## this function is called by permatswap
`permatswap1` <-
function(m, method="quasiswap", fixedmar="both", shuffle="both", strata=NULL,
         mtype="count", thin = 1)
{
    ## internal functions
    indshuffle <- function(x)
    {
       N <- length(x)
       n <- sum(x)
       out <- numeric(N)
       names(out) <- 1:N
       y <- table(sample(1:N, n, replace = TRUE))
       out[names(out) %in% names(y)] <- y
       names(out) <- NULL
       out
    }
    bothshuffle <- function(x, y=1)
    {
        x[x!=0] <- indshuffle(x[x!=0] - y) + y
        sample(x)
    }
    if (!identical(all.equal(m, round(m)), TRUE))
       stop("function accepts only integers (counts)")
    mtype <- match.arg(mtype, c("prab", "count"))
    fixedmar <- match.arg(fixedmar, c("rows", "columns", "both"))
    shuffle <- match.arg(shuffle, c("samp", "both"))
    count <- mtype == "count"
    if (thin < 1)
        stop("'thin' must be >= 1")
    if (count) {
        method <- match.arg(method, c("swap", "quasiswap", "swsh", "abuswap"))
        if (method == "swap") {
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
                if (fixedmar != "both")
                    stop("'fixedmar' must be \"both\"")
            }
        }
    } else {
        if (fixedmar != "both")
            stop("if 'mtype=\"prab\"', 'fixedmar' must be \"both\"")
        method <- match.arg(method, c("swap", "quasiswap", "tswap", "backtracking"))
        isSeq <- method != "quasiswap"
        }

    m <- as.matrix(m)
    att <- attributes(m)
    n.row <- nrow(m)
    n.col <- ncol(m)
    if (mtype == "prab") 
        m <- ifelse(m > 0, 1, 0)

    if (is.null(strata))
        str <- as.factor(rep(1, n.row))
        else str <- as.factor(strata)[drop = TRUE]

    levels(str) <- 1:length(unique(str))
    str <- as.numeric(str)
    nstr <- length(unique(str))
    if (any(tapply(str,list(str),length) == 1))
        stop("strata should contain at least 2 observations")

    perm <- matrix(0, n.row, n.col)

    for (j in 1:nstr) {
        id <- which(str == j)
        temp <- m[id,]
        nn.row <- nrow(m[id,])
        nn.col <- ncol(m[id,])
        if (isSeq) {
            if (count) {
                if (method == "swap")
                    perm[id,] <- .C("swapcount",
                                m = as.double(temp),
                                as.integer(nn.row),
                                as.integer(nn.col),
                                as.integer(thin),
                                PACKAGE = "vegan")$m
                if (method == "abuswap")
                    perm[id,] <- .C("abuswap",
                                m = as.double(temp),
                                as.integer(nn.row),
                                as.integer(nn.col),
                                as.integer(thin),
                                as.integer(direct),
                                PACKAGE = "vegan")$m
            } else {
                perm[id,] <- commsimulator(temp, method=method, thin=thin)
            }
        } else {
            if (method != "swsh") {
                r2tabs <- r2dtable(1, rowSums(m[id,]), colSums(m[id,]))[[1]]
            } else {
                tempPos <- temp[temp > 0]
            }
            if (count) {
                if (method != "swsh") {
                    ms <- sum(m[id,] > 0)
                    tmp <- r2tabs
                    ## if fills are equal, no need to restore fill
                    if (sum(tmp > 0) != ms) {
                        tmp <- .C("rswapcount",
                                    m = as.double(tmp),
                                    as.integer(nn.row),
                                    as.integer(nn.col),
                                    as.integer(ms),
                                    PACKAGE="vegan")$m
                    }
                    perm[id,] <- matrix(tmp, nrow(perm[id,]), ncol(perm[id,]))
                } else { # method == "swsh"
                    tmp <- commsimulator(temp, method="quasiswap")
                    if (shuffle == "samp") {
                        tmp[tmp > 0] <- sample(tempPos)
                    } else {
                        tmp[tmp > 0] <- bothshuffle(tempPos)
                    }
                    perm[id,] <- tmp
                }
            } else perm[id,] <- commsimulator(temp, method=method)
        }
    } # for j end
    perm
}


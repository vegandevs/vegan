## this lists all known algos in vegan and more
## if method is commsim object, it is returned
## if it is character, switch returns the right one, else stop with error
## so it can be used instead of match.arg(method) in other functions
## NOTE: very very long -- but it can be a central repository of algos
## NOTE 2: storage mode coercions are avoided here
## (with no apparent effect on speed), it should be
## handled by nullmodel and commsim characteristics
make.commsim <-
function(method)
{
    algos <- list(
        "r00" = commsim(method="r00", binary=TRUE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- matrix(0L, nr * nc, n)
            for (k in seq_len(n))
                out[sample.int(nr * nc, s), k] <- 1L
            dim(out) <- c(nr, nc, n)
            out
        }),
        "c0" = commsim(method="c0", binary=TRUE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- array(0L, c(nr, nc, n))
            J <- seq_len(nc)
            for (k in seq_len(n))
                for (j in J)
                    out[sample.int(nr, cs[j]), j, k] <- 1L
            out
        }),
        "r0" = commsim(method="r0", binary=TRUE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- array(0L, c(nr, nc, n))
            I <- seq_len(nr)
            for (k in seq_len(n))
                for (i in I)
                    out[i, sample.int(nc, rs[i]), k] <- 1L
            out
        }),
        "r0_old" = commsim(method="r0_old", binary=TRUE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- array(0L, c(nr, nc, n))
            I <- seq_len(nr)
            p <- rep(1, nc)
            for (k in seq_len(n))
                for (i in I)
                    out[i, sample.int(nc, rs[i], prob = p), k] <- 1L
            out
        }),
        "r1" = commsim(method="r1", binary=TRUE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- array(0L, c(nr, nc, n))
            I <- seq_len(nr)
            storage.mode(cs) <- "double"
            for (k in seq_len(n))
                for (i in I)
                    out[i, sample.int(nc, rs[i], prob=cs), k] <- 1L
            out
        }),
        "r2" = commsim(method="r2", binary=TRUE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- array(0L, c(nr, nc, n))
            p <- cs * cs
            I <- seq_len(nr)
            for (k in seq_len(n))
                for (i in I)
                    out[i, sample.int(nc, rs[i], prob=p), k] <- 1L
            out
        }),
        "quasiswap" = commsim(method="quasiswap", binary=TRUE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- array(unlist(r2dtable(n, rs, cs)), c(nr, nc, n))
            storage.mode(out) <- "integer"
            for (k in seq_len(n))
                out[,,k] <- .C("quasiswap",
                    m = out[,,k], nr, nc, PACKAGE = "vegan")$m
            out
        }),
        "swap" = commsim(method="swap", binary=TRUE, isSeq=TRUE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- array(0L, c(nr, nc, n))
            out[,,1] <- .C("swap",
                m = x, nr, nc, thin, PACKAGE = "vegan")$m
            for (k in seq_len(n-1))
                out[,,k+1] <- .C("swap",
                    m = out[,,k], nr, nc, thin,
                    PACKAGE = "vegan")$m
            out
        }),
        "tswap" = commsim(method="tswap", binary=TRUE, isSeq=TRUE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- array(0L, c(nr, nc, n))
            out[,,1] <- .C("trialswap",
                m = x, nr, nc, thin, PACKAGE = "vegan")$m
            for (k in seq_len(n-1))
                out[,,k+1] <- .C("trialswap",
                    m = out[,,k], nr, nc, thin, PACKAGE = "vegan")$m
            out
        }),
        "backtrack" = commsim(method="backtrack", binary=TRUE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            btrfun <- function() {
                all <- matrix(as.integer(1:(nr * nc)), nrow = nr, ncol = nc)
                out <- matrix(0L, nrow = nr, ncol = nc)
                free <- matrix(as.integer(1:(nr * nc)), nrow = nr)
                icount <- integer(length(rs))
                jcount <- integer(length(cs))
                prob <- outer(rs, cs, "*")
                ij <- sample(free, prob = prob)
                i <- (ij - 1)%%nr + 1
                j <- (ij - 1)%/%nr + 1
                for (k in seq_along(ij)) {
                    if (icount[i[k]] < rs[i[k]] && jcount[j[k]] < cs[j[k]]) {
                        out[ij[k]] <- 1L
                        icount[i[k]] <- icount[i[k]] + 1L
                        jcount[j[k]] <- jcount[j[k]] + 1L
                    }
                }
                ndrop <- 1
                for (i in seq_len(10000)) {
                    oldout <- out
                    oldn <- sum(out)
                    drop <- sample(all[out == 1L], ndrop)
                    out[drop] <- 0L
                    candi <- outer(rowSums(out) < rs, colSums(out) < cs, "&") & out == 0L
                    while (sum(candi) > 0) {
                        if (sum(candi) > 1)
                          ij <- sample(all[candi], 1)
                        else ij <- all[candi]
                        out[ij] <- 1L
                        candi <- outer(rowSums(out) < rs, colSums(out) < cs, "&") & out == 0
                    }
                    if (sum(out) >= fill)
                        break
                    if (oldn >= sum(out))
                        ndrop <- min(ndrop + 1, 4)
                    else ndrop <- 1
                    if (oldn > sum(out))
                        out <- oldout
                }
                out
            }
            out <- array(0L, c(nr, nc, n))
            for (k in seq_len(n))
                out[, , k] <- btrfun()
            out
        }),
        "r2dtable" = commsim(method="r2dtable", binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            out <- array(unlist(r2dtable(n, rs, cs)), c(nr, nc, n))
            storage.mode(out) <- "integer"
            out
        }),
        "swap_count" = commsim(method="swap_count", binary=FALSE, isSeq=TRUE,
        mode="integer",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            out <- array(0L, c(nr, nc, n))
            out[,,1] <- .C("swapcount",
                m = x, nr, nc, thin, PACKAGE = "vegan")$m
            for (k in seq_len(n-1))
                out[,,k+1] <- .C("swapcount",
                    m = out[,,k], nr, nc, thin, PACKAGE = "vegan")$m
            out
        }),
        "quasiswap_count" = commsim(method="quasiswap_count", binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            out <- array(unlist(r2dtable(n, rs, cs)), c(nr, nc, n))
            storage.mode(out) <- "integer"
            for (k in seq_len(n))
                out[,,k] <- .C("rswapcount",
                    m = out[,,k], nr, nc, fill, PACKAGE = "vegan")$m
            out
        }),
        "swsh_samp" = commsim(method="swsh_samp", binary=FALSE, isSeq=FALSE,
        mode="double",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            nz <- x[x > 0]
            out <- array(unlist(r2dtable(fill, rf, cf)), c(nr, nc, n))
            storage.mode(out) <- "double"
            for (k in seq_len(n)) {
                out[,,k] <- .C("quasiswap",
                    m = as.integer(out[,,k]), nr, nc, PACKAGE = "vegan")$m
                out[,,k][out[,,k] > 0] <- sample(nz) # we assume that length(nz)>1
            }
            out
        }),
        "swsh_both" = commsim(method="swsh_both", binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            indshuffle <- function(x) {
                drop(rmultinom(1, sum(x), rep(1, length(x))))
            }
            nz <- as.integer(x[x > 0])
            out <- array(unlist(r2dtable(fill, rf, cf)), c(nr, nc, n))
            storage.mode(out) <- "integer"
            for (k in seq_len(n)) {
                out[,,k] <- .C("quasiswap",
                    m = out[,,k], nr, nc, PACKAGE = "vegan")$m
                out[,,k][out[,,k] > 0] <- indshuffle(nz - 1L) + 1L  # we assume that length(nz)>1
            }
            out
        }),
        "swsh_samp_r" = commsim(method="swsh_samp_r", binary=FALSE, isSeq=FALSE,
        mode="double",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            out <- array(unlist(r2dtable(fill, rf, cf)), c(nr, nc, n))
            storage.mode(out) <- "double"
            I <- seq_len(nr)
            for (k in seq_len(n)) {
                out[,,k] <- .C("quasiswap",
                    m = as.integer(out[,,k]), nr, nc, PACKAGE = "vegan")$m
                for (i in I) {
                    nz <- x[i,][x[i,] > 0]
                    if (length(nz) == 1)
                        out[i,,k][out[i,,k] > 0] <- nz
                    if (length(nz) > 1)
                        out[i,,k][out[i,,k] > 0] <- sample(nz)
                }
            }
            out
        }),
        "swsh_samp_c" = commsim(method="swsh_samp_c", binary=FALSE, isSeq=FALSE,
        mode="double",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            out <- array(unlist(r2dtable(fill, rf, cf)), c(nr, nc, n))
            storage.mode(out) <- "double"
            J <- seq_len(nc)
            for (k in seq_len(n)) {
                out[,,k] <- .C("quasiswap",
                    m = as.integer(out[,,k]), nr, nc, PACKAGE = "vegan")$m
                for (j in J) {
                    nz <- x[,j][x[,j] > 0]
                    if (length(nz) == 1)
                        out[,j,k][out[,j,k] > 0] <- nz
                    if (length(nz) > 1)
                        out[,j,k][out[,j,k] > 0] <- sample(nz)
                }
            }
            out
        }),
        "swsh_both_r" = commsim(method="swsh_both_r", binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            indshuffle <- function(x) {
                drop(rmultinom(1, sum(x), rep(1, length(x))))
            }
            I <- seq_len(nr)
            out <- array(unlist(r2dtable(fill, rf, cf)), c(nr, nc, n))
            storage.mode(out) <- "integer"
            for (k in seq_len(n)) {
                out[,,k] <- .C("quasiswap",
                    m = out[,,k], nr, nc, PACKAGE = "vegan")$m
                for (i in I) {
                    nz <- as.integer(x[i,][x[i,] > 0])
                    if (length(nz) == 1)
                        out[i,,k][out[i,,k] > 0] <- nz
                    if (length(nz) > 1)
                        out[i,,k][out[i,,k] > 0] <- indshuffle(nz - 1L) + 1L
                }
            }
            out
        }),
        "swsh_both_c" = commsim(method="swsh_both_c", binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            indshuffle <- function(x) {
                drop(rmultinom(1, sum(x), rep(1, length(x))))
            }
            J <- seq_len(nc)
            out <- array(unlist(r2dtable(fill, rf, cf)), c(nr, nc, n))
            storage.mode(out) <- "integer"
            for (k in seq_len(n)) {
                out[,,k] <- .C("quasiswap",
                    m = out[,,k], nr, nc,  PACKAGE = "vegan")$m
                for (j in J) {
                    nz <- as.integer(x[,j][x[,j] > 0])
                    if (length(nz) == 1)
                        out[,j,k][out[,j,k] > 0] <- nz
                    if (length(nz) > 1)
                        out[,j,k][out[,j,k] > 0] <- indshuffle(nz - 1L) + 1L
                }
            }
            out
        }),
        "abuswap_r" = commsim(method="abuswap_r", binary=FALSE, isSeq=TRUE,
        mode="double",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            out <- array(0, c(nr, nc, n))
            out[,,1] <- .C("abuswap",
                m = x, nr, nc, thin, 1L, PACKAGE = "vegan")$m
            for (k in seq_len(n-1))
                out[,,k+1] <- .C("abuswap",
                    m = out[,,k], nr, nc, thin, 1L, PACKAGE = "vegan")$m
            out
        }),
        "abuswap_c" = commsim(method="abuswap_c", binary=FALSE, isSeq=TRUE,
        mode="double",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            out <- array(0, c(nr, nc, n))
            out[,,1] <- .C("abuswap",
                m = x, nr, nc, thin, 0L, PACKAGE = "vegan")$m
            for (k in seq_len(n-1))
                out[,,k+1] <- .C("abuswap",
                    m = out[,,k], nr, nc, thin, 0L, PACKAGE = "vegan")$m
            out
        }),
        "r00_samp" = commsim(method="r00_samp", binary=FALSE, isSeq=FALSE,
        mode="double",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- matrix(0, nr * nc, n)
            for (k in seq_len(n))
                out[, k] <- sample(x)
            dim(out) <- c(nr, nc, n)
            out
        }),
        "c0_samp" = commsim(method="c0_samp", binary=FALSE, isSeq=FALSE,
        mode="double",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- array(0, c(nr, nc, n))
            J <- seq_len(nc)
            for (k in seq_len(n))
                for (j in J)
                    out[, j, k] <- sample(x[,j])
            out
        }),
        "r0_samp" = commsim(method="r0_samp", binary=FALSE, isSeq=FALSE,
        mode="double",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- array(0, c(nr, nc, n))
            I <- seq_len(nr)
            for (k in seq_len(n))
                for (i in I)
                    out[i, , k] <- sample(x[i,])
            out
        }),
        "r00_ind" = commsim(method="r00_ind", binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            indshuffle <- function(x) {
                drop(rmultinom(1, sum(x), rep(1, length(x))))
            }
            out <- matrix(0L, nr * nc, n)
            for (k in seq_len(n))
                out[, k] <- indshuffle(x)
            dim(out) <- c(nr, nc, n)
            out
        }),
        "c0_ind" = commsim(method="c0_ind", binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            indshuffle <- function(x) {
                drop(rmultinom(1, sum(x), rep(1, length(x))))
            }
            out <- array(0L, c(nr, nc, n))
            J <- seq_len(nc)
            for (k in seq_len(n))
                for (j in J)
                    out[, j, k] <- indshuffle(x[,j])
            out
        }),
        "r0_ind" = commsim(method="r0_ind", binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            indshuffle <- function(x) {
                drop(rmultinom(1, sum(x), rep(1, length(x))))
            }
            out <- array(0L, c(nr, nc, n))
            I <- seq_len(nr)
            for (k in seq_len(n))
                for (i in I)
                    out[i, , k] <- indshuffle(x[i,])
            out
        }),
        "r00_both" = commsim(method="r00_both", binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            indshuffle <- function(x) {
                drop(rmultinom(1, sum(x), rep(1, length(x))))
            }
            out <- matrix(0L, nr * nc, n)
            for (k in seq_len(n)) {
                out[,k][x > 0] <- indshuffle(x[x > 0] - 1L) + 1L
                out[,k] <- sample(out[,k])
            }
            dim(out) <- c(nr, nc, n)
            out
        }),
        "c0_both" = commsim(method="c0_both", binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            indshuffle <- function(x) {
                drop(rmultinom(1, sum(x), rep(1, length(x))))
            }
            out <- array(0L, c(nr, nc, n))
            J <- seq_len(nc)
            for (k in seq_len(n))
                for (j in J) {
                    out[,j,k][x[,j] > 0] <- indshuffle(x[,j][x[,j] > 0] - 1L) + 1L
                    out[,j,k] <- sample(out[,j,k])
                }
            out
        }),
        "r0_both" = commsim(method="r0_both", binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            indshuffle <- function(x) {
                drop(rmultinom(1, sum(x), rep(1, length(x))))
            }
            out <- array(0L, c(nr, nc, n))
            I <- seq_len(nr)
            for (k in seq_len(n))
                for (i in I) {
                    out[i,,k][x[i,] > 0] <- indshuffle(x[i,][x[i,] > 0] - 1L) + 1L
                    out[i,,k] <- sample(out[i,,k])
                }
            out
        })
    )
    if (missing(method))
        return(names(algos))
    if (inherits(method, "commsim"))
        return(method)
    method <- match.arg(method, sort(names(algos)))
    algos[[method]]
}

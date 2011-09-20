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
    if (inherits(method, "commsim"))
        return(method)
    switch(method, 
        "r00" = return(commsim(method=method, binary=TRUE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- matrix(0L, nr * nc, n)
            for (k in seq_len(n))
                out[sample.int(nr * nc, s), k] <- 1
            dim(out) <- c(nr, nc, n)
            out
        })),
        "c0" = return(commsim(method=method, binary=TRUE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- array(0L, c(nr, nc, n))
            J <- seq_len(nc)
            for (k in seq_len(n))
                for (j in J)
                    out[sample.int(nr, cs[j]), j, k] <- 1
            out
        })),
        "r0" = return(commsim(method=method, binary=TRUE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- array(0L, c(nr, nc, n))
            I <- seq_len(nr)
            for (k in seq_len(n))
                for (i in I)
                    out[i, sample.int(nc, rs[i]), k] <- 1
            out
        })),
        "r1" = return(commsim(method=method, binary=TRUE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- array(0L, c(nr, nc, n))
            I <- seq_len(nr)
            for (k in seq_len(n))
                for (i in I)
                    out[i, sample.int(nc, rs[i], prob=cs), k] <- 1
            out
        })),
        "r2" = return(commsim(method=method, binary=TRUE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- array(0L, c(nr, nc, n))
            p <- cs * cs
            I <- seq_len(nr)
            for (k in seq_len(n))
                for (i in I)
                    out[i, sample.int(nc, rs[i], prob=p), k] <- 1
            out
        })),
        "quasiswap" = return(commsim(method=method, binary=TRUE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- array(unlist(r2dtable(n, rs, cs)), c(nr, nc, n))
            storage.mode(out) <- "integer"
            for (k in seq_len(n))
                out[,,k] <- .C("quasiswap", 
                    m = out[,,k], nr, nc, PACKAGE = "vegan")$m
            out
        })),
        "swap" = return(commsim(method=method, binary=TRUE, isSeq=TRUE,
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
        })),
        "tswap" = return(commsim(method=method, binary=TRUE, isSeq=TRUE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- array(0L, c(nr, nc, n))
            out[,,1] <- .C("trialswap", 
                m = x, nr, nc, thin, PACKAGE = "vegan")$m
            for (k in seq_len(n-1))
                out[,,k+1] <- .C("trialswap", 
                    m = out[,,k], nr, nc, thin, PACKAGE = "vegan")$m
            out
        })),
        "backtrack" = return(commsim(method=method, binary=TRUE, isSeq=FALSE,
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
                for (k in 1:length(ij)) {
                    if (icount[i[k]] < rs[i[k]] && jcount[j[k]] < cs[j[k]]) {
                        out[ij[k]] <- 1
                        icount[i[k]] <- icount[i[k]] + 1
                        jcount[j[k]] <- jcount[j[k]] + 1
                    }
                }
                ndrop <- 1
                for (i in 1:10000) {
                    oldout <- out
                    oldn <- sum(out)
                    drop <- sample(all[out == 1], ndrop)
                    out[drop] <- 0
                    candi <- outer(rowSums(out) < rs, colSums(out) < cs, "&") & out == 0
                    while (sum(candi) > 0) {
                        if (sum(candi) > 1) 
                          ij <- sample(all[candi], 1)
                        else ij <- all[candi]
                        out[ij] <- 1
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
            out <- array(0, c(nr, nc, n))
            for (k in seq_len(n))
                out[, , k] <- btrfun()
            out
        })),
        "r2dtable" = return(commsim(method=method, binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            out <- array(unlist(r2dtable(n, rs, cs)), c(nr, nc, n))
            storage.mode(out) <- "integer"
            out
        })),
        "swap_count" = return(commsim(method=method, binary=FALSE, isSeq=TRUE,
        mode="integer",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            out <- array(0L, c(nr, nc, n))
            out[,,1] <- .C("swapcount", 
#                m = as.double(x), nr, nc, thin, PACKAGE = "vegan")$m
                m = x, nr, nc, thin, PACKAGE = "vegan")$m
            for (k in seq_len(n-1))
                out[,,k+1] <- .C("swapcount", 
#                    m = as.double(out[,,k]), nr, nc, thin, PACKAGE = "vegan")$m
                    m = out[,,k], nr, nc, thin, PACKAGE = "vegan")$m
            out
        })),
        "quasiswap_count" = return(commsim(method=method, binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            out <- array(unlist(r2dtable(n, rs, cs)), c(nr, nc, n))
            storage.mode(out) <- "integer"
            for (k in seq_len(n))
                out[,,k] <- .C("rswapcount", 
#                    m = as.double(out[,,k]), nr, nc, fill, PACKAGE = "vegan")$m
                    m = out[,,k], nr, nc, fill, PACKAGE = "vegan")$m
            out
        })),
        "swsh_samp" = return(commsim(method=method, binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            nz <- as.integer(x[x > 0])
            out <- array(unlist(r2dtable(fill, rf, cf)), c(nr, nc, n))
            storage.mode(out) <- "integer"
            for (k in seq_len(n)) {
                out[,,k] <- .C("quasiswap", 
                    m = out[,,k], nr, nc, PACKAGE = "vegan")$m
                out[,,k][out[,,k] > 0] <- sample(nz)
            }
            out
        })),
        "swsh_both" = return(commsim(method=method, binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            indshuffle <- function(x) {
                out <- integer(length(x))
                tmp <- table(sample.int(length(x), sum(x), replace = TRUE))
                out[as.integer(names(tmp))] <- as.integer(tmp)
                out
            }
            nz <- as.integer(x[x > 0])
            out <- array(unlist(r2dtable(fill, rf, cf)), c(nr, nc, n))
            storage.mode(out) <- "integer"
            for (k in seq_len(n)) {
                out[,,k] <- .C("quasiswap", 
                    m = out[,,k], nr, nc, PACKAGE = "vegan")$m
                out[,,k][out[,,k] > 0] <- sample(indshuffle(nz - 1L) + 1L)
            }
            out
        })),
        "swsh_samp_r" = return(commsim(method=method, binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            out <- array(unlist(r2dtable(fill, rf, cf)), c(nr, nc, n))
            storage.mode(out) <- "integer"
            I <- seq_len(nr)
            for (k in seq_len(n)) {
                out[,,k] <- .C("quasiswap", 
                    m = out[,,k], nr, nc, PACKAGE = "vegan")$m
                for (i in I)
                    out[i,,k][out[i,,k] > 0] <- sample(as.integer(x[i,][x[i,] > 0]))
            }
            out
        })),
        "swsh_samp_c" = return(commsim(method=method, binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            out <- array(unlist(r2dtable(fill, rf, cf)), c(nr, nc, n))
            storage.mode(out) <- "integer"
            J <- seq_len(nc)
            for (k in seq_len(n)) {
                out[,,k] <- .C("quasiswap", 
                    m = out[,,k], nr, nc, PACKAGE = "vegan")$m
                for (j in J)
                    out[,j,k][out[,j,k] > 0] <- sample(as.integer(x[,j][x[,j] > 0]))
            }
            out
        })),
        "swsh_both_r" = return(commsim(method=method, binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            indshuffle <- function(x) {
                out <- integer(length(x))
                tmp <- table(sample.int(length(x), sum(x), replace = TRUE))
                out[as.integer(names(tmp))] <- as.integer(tmp)
                out
            }
            I <- seq_len(nr)
            out <- array(unlist(r2dtable(fill, rf, cf)), c(nr, nc, n))
            storage.mode(out) <- "integer"
            for (k in seq_len(n)) {
                out[,,k] <- .C("quasiswap", 
                    m = out[,,k], nr, nc, PACKAGE = "vegan")$m
                for (i in I)
                    out[i,,k][out[i,,k] > 0] <- sample(indshuffle(as.integer(x[i,][x[i,] > 0]) - 1L) + 1L)
            }
            out
        })),
        "swsh_both_c" = return(commsim(method=method, binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            indshuffle <- function(x) {
                out <- integer(length(x))
                tmp <- table(sample.int(length(x), sum(x), replace = TRUE))
                out[as.integer(names(tmp))] <- as.integer(tmp)
                out
            }
            J <- seq_len(nc)
            out <- array(unlist(r2dtable(fill, rf, cf)), c(nr, nc, n))
            storage.mode(out) <- "integer"
            for (k in seq_len(n)) {
                out[,,k] <- .C("quasiswap", 
                    m = out[,,k], nr, nc,  PACKAGE = "vegan")$m
                for (j in J)
                    out[,j,k][out[,j,k] > 0] <- sample(indshuffle(as.integer(x[,j][x[,j] > 0]) - 1L) + 1L)
            }
            out
        })),
        "abuswap_r" = return(commsim(method=method, binary=FALSE, isSeq=TRUE,
        mode="double",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            out <- array(0, c(nr, nc, n))
            out[,,1] <- .C("abuswap", 
                m = x, nr, nc, thin, 1L, PACKAGE = "vegan")$m
            for (k in seq_len(n-1))
                out[,,k+1] <- .C("abuswap", 
                    m = out[,,k], nr, nc, thin, 1L, PACKAGE = "vegan")$m
            out
        })),
        "abuswap_c" = return(commsim(method=method, binary=FALSE, isSeq=TRUE,
        mode="double",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            out <- array(0, c(nr, nc, n))
            out[,,1] <- .C("abuswap", 
                m = x, nr, nc, thin, 0L, PACKAGE = "vegan")$m
            for (k in seq_len(n-1))
                out[,,k+1] <- .C("abuswap", 
                    m = out[,,k], nr, nc, thin, 0L, PACKAGE = "vegan")$m
            out
        })),
        "r00_samp" = return(commsim(method=method, binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- matrix(0L, nr * nc, n)
            for (k in seq_len(n))
                out[, k] <- sample(x)
            dim(out) <- c(nr, nc, n)
            out
        })),
        "c0_samp" = return(commsim(method=method, binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- array(0L, c(nr, nc, n))
            J <- seq_len(nc)
            for (k in seq_len(n))
                for (j in J)
                    out[, j, k] <- sample(x[,j])
            out
        })),
        "r0_samp" = return(commsim(method=method, binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- array(0L, c(nr, nc, n))
            I <- seq_len(nr)
            for (k in seq_len(n))
                for (i in I)
                    out[i, , k] <- sample(x[i,])
            out
        })),
        "r00_ind" = return(commsim(method=method, binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            indshuffle <- function(x) {
                out <- integer(length(x))
                tmp <- table(sample.int(length(x), sum(x), replace = TRUE))
                out[as.integer(names(tmp))] <- as.integer(tmp)
                out
            }
            out <- matrix(0L, nr * nc, n)
            for (k in seq_len(n))
                out[, k] <- indshuffle(x)
            dim(out) <- c(nr, nc, n)
            out
        })),
        "c0_ind" = return(commsim(method=method, binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            indshuffle <- function(x) {
                out <- integer(length(x))
                tmp <- table(sample.int(length(x), sum(x), replace = TRUE))
                out[as.integer(names(tmp))] <- as.integer(tmp)
                out
            }
            out <- array(0L, c(nr, nc, n))
            J <- seq_len(nc)
            for (k in seq_len(n))
                for (j in J)
                    out[, j, k] <- indshuffle(x[,j])
            out
        })),
        "r0_ind" = return(commsim(method=method, binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            indshuffle <- function(x) {
                out <- integer(length(x))
                tmp <- table(sample.int(length(x), sum(x), replace = TRUE))
                out[as.integer(names(tmp))] <- as.integer(tmp)
                out
            }
            out <- array(0L, c(nr, nc, n))
            I <- seq_len(nr)
            for (k in seq_len(n))
                for (i in I)
                    out[i, , k] <- indshuffle(x[i,])
            out
        })),
        "r00_both" = return(commsim(method=method, binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            indshuffle <- function(x) {
                out <- integer(length(x))
                tmp <- table(sample.int(length(x), sum(x), replace = TRUE))
                out[as.integer(names(tmp))] <- as.integer(tmp)
                out
            }
            out <- matrix(0L, nr * nc, n)
            for (k in seq_len(n)) {
                out[,k][x > 0] <- indshuffle(x[x > 0] - 1L) + 1L
                out[,k] <- sample(out[,k])
            }
            dim(out) <- c(nr, nc, n)
            out
        })),
        "c0_both" = return(commsim(method=method, binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            indshuffle <- function(x) {
                out <- integer(length(x))
                tmp <- table(sample.int(length(x), sum(x), replace = TRUE))
                out[as.integer(names(tmp))] <- as.integer(tmp)
                out
            }
            out <- array(0L, c(nr, nc, n))
            J <- seq_len(nc)
            for (k in seq_len(n))
                for (j in J) {
                    out[,j,k][x[,j] > 0] <- indshuffle(x[,j][x[,j] > 0] - 1L) + 1L
                    out[,j,k] <- sample(out[,j,k])
                }
            out
        })),
        "r0_both" = return(commsim(method=method, binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            indshuffle <- function(x) {
                out <- integer(length(x))
                tmp <- table(sample.int(length(x), sum(x), replace = TRUE))
                out[as.integer(names(tmp))] <- as.integer(tmp)
                out
            }
            out <- array(0L, c(nr, nc, n))
            I <- seq_len(nr)
            for (k in seq_len(n))
                for (i in I) {
                    out[i,,k][x[i,] > 0] <- indshuffle(x[i,][x[i,] > 0] - 1L) + 1L
                    out[i,,k] <- sample(out[i,,k])
                }
            out
        }))
    )
    stop("\"", method, "\" method not found")
}

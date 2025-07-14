## this lists all known algos in vegan and more
## if method is commsim object, it is returned
## if it is character, switch returns the right one, else stop with error
## so it can be used instead of match.arg(method) in other functions
## NOTE: very very long -- but it can be a central repository of algos
## NOTE 2: storage mode coercions are avoided here
## (with no apparent effect on speed), it should be
## handled by nullmodel and commsim characteristics
`make.commsim` <-
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
            if (nr < 2L || nc < 2L)
                stop(sQuote(method)," needs at least two items", call. = FALSE)
            out <- array(unlist(r2dtable(n, rs, cs)), c(nr, nc, n))
            storage.mode(out) <- "integer"
            .Call(do_qswap, out, n, thin, "quasiswap", PACKAGE = "vegan")
        }),
        "greedyqswap" = commsim(method="greedyqswap", binary=TRUE,
               isSeq=FALSE, mode="integer",
               fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            if (nr < 2L || nc < 2L)
                stop(sQuote(method)," needs at least two items", call. = FALSE)
            out <- array(unlist(r2dtable(n, rs, cs)), c(nr, nc, n))
            storage.mode(out) <- "integer"
            .Call(do_greedyqswap, out, n, thin, fill, PACKAGE = "vegan")
        }),
        "swap" = commsim(method="swap", binary = TRUE, isSeq=TRUE,
        mode = "integer",
        fun = function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            if (nr < 2L || nc < 2L)
                stop(sQuote(method)," needs at least two items", call. = FALSE)
            ## no checkerboard 2x2 matrices: infinite loop
            if (nestedchecker(x)$statistic == 0)
                stop(sQuote(method), " needs checkerboard data: check with nestedchecker(x)", call. = FALSE)
            .Call(do_swap, as.matrix(x), n, thin, "swap", PACKAGE = "vegan")
        }),
        "tswap" = commsim(method="tswap", binary = TRUE, isSeq=TRUE,
        mode = "integer",
        fun = function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            if (nr < 2L || nc < 2L)
                stop(sQuote(method)," needs at least two items", call. = FALSE)
            .Call(do_swap, as.matrix(x), n, thin, "trialswap", PACKAGE = "vegan")
        }),
        "curveball" = commsim(method="curveball", binary = TRUE, isSeq=TRUE,
        mode = "integer",
        fun = function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            if (nrow(x) < 2)
                stop(sQuote(method), " needs at least two rows", call. = FALSE)
            .Call(do_curveball, as.matrix(x), n, thin, PACKAGE = "vegan")
        }),
        "backtrack" = commsim(method="backtrack", binary = TRUE,
                               isSeq = FALSE, mode = "integer",
        fun = function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            .Call(do_backtrack, n, rs, cs, PACKAGE = "vegan")
        }),
        "r2dtable" = commsim(method="r2dtable", binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            if (nr < 2L || nc < 2L)
                stop(sQuote(method)," needs at least two items", call. = FALSE)
            out <- array(unlist(r2dtable(n, rs, cs)), c(nr, nc, n))
            storage.mode(out) <- "integer"
            out
        }),
        "swap_count" = commsim(method="swap_count", binary = FALSE,
        isSeq=TRUE, mode = "integer",
        fun = function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            if (nr < 2L || nc < 2L)
                stop(sQuote(method)," needs at least two items", call. = FALSE)
            ## no checkerboard 2x2 matrices: infinite loop
            if (nestedchecker(x)$statistic == 0)
                stop(sQuote(method), " needs checkerboard data: check with nestedchecker(x)", call. = FALSE)
            .Call(do_swap, as.matrix(x), n, thin, "swapcount", PACKAGE = "vegan")
        }),
        "quasiswap_count" = commsim(method="quasiswap_count", binary=FALSE,
        isSeq=FALSE, mode="integer",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            if (nr < 2L || nc < 2L)
                stop(sQuote(method)," needs at least two items", call. = FALSE)
            ## no checkerboard 2x2 matrices: infinite loop
            if (nestedchecker(x)$statistic == 0)
                stop(sQuote(method), " needs checkerboard data: check with nestedchecker(x)", call. = FALSE)
            out <- array(unlist(r2dtable(n, rs, cs)), c(nr, nc, n))
            storage.mode(out) <- "integer"
            .Call(do_qswap, out, n, fill, "rswapcount", PACKAGE = "vegan")
        }),
        "swsh_samp" = commsim(method="swsh_samp", binary=FALSE, isSeq=FALSE,
        mode="double",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            if (nr < 2L || nc < 2L)
                stop(sQuote(method)," needs at least two items", call. = FALSE)
            nz <- x[x > 0]
            out <- array(unlist(r2dtable(n, rf, cf)), c(nr, nc, n))
            ## do_qswap changes 'out' within the function
            .Call(do_qswap, out, n, thin, "quasiswap", PACKAGE = "vegan")
            storage.mode(out) <- "double"
            for (k in seq_len(n)) {
                out[,,k][out[,,k] > 0] <- sample(nz) # we assume that length(nz)>1
            }
            out
        }),
        "swsh_both" = commsim(method="swsh_both", binary=FALSE, isSeq=FALSE,
        mode="integer",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            if (nr < 2L || nc < 2L)
                stop(sQuote(method)," needs at least two items", call. = FALSE)
            indshuffle <- function(x) {
                drop(rmultinom(1, sum(x), rep(1, length(x))))
            }
            nz <- as.integer(x[x > 0])
            out <- array(unlist(r2dtable(n, rf, cf)), c(nr, nc, n))
            .Call(do_qswap, out, n, thin, "quasiswap", PACKAGE = "vegan")
            storage.mode(out) <- "integer"
            for (k in seq_len(n)) {
                out[,,k][out[,,k] > 0] <- indshuffle(nz - 1L) + 1L  # we assume that length(nz)>1
            }
            out
        }),
        "swsh_samp_r" = commsim(method="swsh_samp_r", binary=FALSE, isSeq=FALSE,
        mode="double",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            if (nr < 2L || nc < 2L)
                stop(sQuote(method)," needs at least two items", call. = FALSE)
            out <- array(unlist(r2dtable(n, rf, cf)), c(nr, nc, n))
            .Call(do_qswap, out, n, thin, "quasiswap", PACKAGE = "vegan")
            storage.mode(out) <- "double"
            I <- seq_len(nr)
            for (k in seq_len(n)) {
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
            if (nr < 2L || nc < 2L)
                stop(sQuote(method)," needs at least two items", call. = FALSE)
            out <- array(unlist(r2dtable(n, rf, cf)), c(nr, nc, n))
            .Call(do_qswap, out, n, thin, "quasiswap", PACKAGE = "vegan")
            storage.mode(out) <- "double"
            J <- seq_len(nc)
            for (k in seq_len(n)) {
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
            if (nr < 2L || nc < 2L)
                stop(sQuote(method)," needs at least two items", call. = FALSE)
            indshuffle <- function(x) {
                drop(rmultinom(1, sum(x), rep(1, length(x))))
            }
            I <- seq_len(nr)
            out <- array(unlist(r2dtable(n, rf, cf)), c(nr, nc, n))
            .Call(do_qswap, out, n, thin, "quasiswap", PACKAGE = "vegan")
            storage.mode(out) <- "integer"
            for (k in seq_len(n)) {
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
            if (nr < 2L || nc < 2L)
                stop(sQuote(method)," needs at least two items", call. = FALSE)
            indshuffle <- function(x) {
                drop(rmultinom(1, sum(x), rep(1, length(x))))
            }
            J <- seq_len(nc)
            out <- array(unlist(r2dtable(n, rf, cf)), c(nr, nc, n))
            .Call(do_qswap, out, n, thin, "quasiswap", PACKAGE = "vegan")
            storage.mode(out) <- "integer"
            for (k in seq_len(n)) {
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
            ## no checkerboard 2x2 matrices: infinite loop
            if (nestedchecker(x)$statistic == 0)
                stop(sQuote(method), " needs checkerboard data: check with nestedchecker(x)", call. = FALSE)
            .Call(do_abuswap, as.matrix(x), n, thin, 1L, PACKAGE = "vegan")
        }),
        "abuswap_c" = commsim(method="abuswap_c", binary=FALSE, isSeq=TRUE,
        mode="double",
        fun=function(x, n, nr, nc, cs, rs, rf, cf, s, fill, thin) {
            ## no checkerboard 2x2 matrices: infinite loop
            if (nestedchecker(x)$statistic == 0)
                stop(sQuote(method), " needs checkerboard data: check with nestedchecker(x)", call. = FALSE)
            .Call(do_abuswap, as.matrix(x), n, thin, 0L, PACKAGE = "vegan")
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
                    out[, j, k] <- if (nr < 2)
                        x[,j] else sample(x[,j])
            out
        }),
        "r0_samp" = commsim(method="r0_samp", binary=FALSE, isSeq=FALSE,
        mode="double",
        fun=function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin) {
            out <- array(0, c(nr, nc, n))
            I <- seq_len(nr)
            for (k in seq_len(n))
                for (i in I)
                    out[i, , k] <- if (nc < 2)
                        x[i,] else sample(x[i,])
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
                    if (sum(x[,j]) > 0) {
                        out[,j,k][x[,j] > 0] <- indshuffle(x[,j][x[,j] > 0] - 1L) + 1L
                        out[,j,k] <- if (nr < 2)
                            out[,j,k] else sample(out[,j,k])
                    }
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
                    if (sum(x[i,]) > 0) {
                        out[i,,k][x[i,] > 0] <- indshuffle(x[i,][x[i,] > 0] - 1L) + 1L
                        out[i,,k] <- if (nc < 2)
                            out[i,,k] else sample(out[i,,k])
                    }
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

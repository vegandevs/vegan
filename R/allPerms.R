`allPerms` <- function(n, control = permControl(), max = 9999,
                       observed = FALSE) {
    ## in-line functions
    `all.free` <- function(n, v = 1:n) {
	if( n == 1 ) {
            matrix(v, 1, 1)
        } else {
            X <- NULL
            for(i in 1:n)
                X <- rbind(X, cbind(v[i],
                                    Recall(n-1, v[-i])))
            X
        }
    }
    `all.series` <- function(n, control) {
        v <- seq_len(n)
        nperms <- numPerms(v, control = control)
	X <- matrix(nrow = nperms, ncol = n)
	for(i in v) {
            X[i,] <- seq(i, length = n)%%n + 1
	}
        ## if mirroring, rev the cols of X[v,]
        ## but only if n > 2
        if(control$mirror && (nperms > 2))
            X[(n+1):(2*n),] <- X[v, rev(v)]
	X
    }
    `all.grid` <- function(n, control) {
        v <- seq_len(n)
        nperms <- numPerms(v, control)
        nr <- control$nrow
        nc <- control$ncol
	X <- matrix(nrow = nperms, ncol = n)
        idx <- 1
        ## ncol == 2 is special case
        if(control$ncol == 2) {
            X <- all.series(n, permControl(type = "series",
                                           mirror = control$mirror,
                                           constant = control$constant)
                            )
        } else {
            for(i in seq_len(nr)) {
                for(j in seq_len(nc)) {
                    ir <- seq(i, length = nr)%%nr
                    ic <- seq(j, length = nc)%%nc
                    ## block 1 - no reversals
                    X[idx, ] <- rep(ic, each = nr) * nr +
                        rep(ir, len = nr * nc) + 1
                    if(control$mirror) {
                        ## block 2 - rev rows but not columns
                        X[idx + n, ] <- rep(ic, each = nr) * nr +
                            rep(rev(ir), len = nr * nc) + 1
                        ## block 3 - rev columns but not rows
                        X[idx + (2*n), ] <- rep(rev(ic), each = nr) *
                            nr + rep(ir, len = nr * nc) + 1
                    }
                    idx <- idx + 1
                }
            }
            if(control$mirror) {
                ## rev columns and rows
                ## no calculations, just rev cols of block 1
                v <- seq_len(n)
                X[((3*n)+1):(4*n), ] <- X[v, rev(v)]
            }
        }
        X
    }
    `all.strata` <- function(n, control) {
        v <- seq_len(n)
        nperms <- numPerms(v, control)
        lev <- length(levels(control$strata))
        X <- matrix(nrow = nperms, ncol = length(control$strata))
        perms <- if(control$type == "free") {
            all.free(lev)
        } else if(control$type == "series") {
            all.series(lev, control = control)
        } else {
            all.grid(lev, control = control)
        }
        sp <- split(v, control$strata)
        for(i in seq_len(nrow(perms)))
            X[i,] <- unname(do.call(c, sp[perms[i,]]))
        X
    }
    ## replacement for recursive function above
    bar <- function(mat, n) {
        res <- vector(mode = "list", length = n)
        for(i in seq_len(n))
            res[[i]] <- mat
        do.call(rbind, res)
    }
    ## start
    v <- n
    ## expand n if a numeric or integer vector of length 1
    if((is.numeric(n) || is.integer(n)) && (length(n) == 1))
         v <- seq_len(n)
    ## number of observations in data
    n <- getNumObs(v)
    ## check permutation scheme and update control
    pcheck <- permCheck(v, control = control, make.all = FALSE)
    control <- pcheck$control
    ## get max number of permutations
    ## originally had:
    ##nperms <- numPerms(v, control = control)
    ## but pcheck contains 'n', the result of call to numPerms
    nperms <- pcheck$n
    ## sanity check - don't let this run away to infinity
    ## esp with type = "free"
    if(nperms > max)
        stop("Number of possible permutations too big (> 'max')")
    type <- control$type
    ##if(type != "strata" && !is.null(control$strata)) {
    if(!control$permute.strata && !is.null(control$strata)) {
        ## permuting within blocks
        ## FIXME: allperms expects samples to be arranged
        ## in order of fac, i.e. all level 1, followed by
        ## all level 2 - fix to allow them to be in any order:
        ## see permuted.index2 for how to do this
        if(control$constant) {
            ## same permutation in each block
            #v <- seq_len(n)
            pg <- unique(table(control$strata))
            control.wi <- permControl(type = control$type,
                                      mirror = control$mirror,
                                      nrow = control$nrow,
                                      ncol = control$ncol)
            nperms <- numPerms(v, control)
            ord <- switch(control$type,
                          free = all.free(pg),
                          series = all.series(pg, control = control.wi),
                          grid = all.grid(pg, control = control.wi))
            perm.wi <- nrow(ord)
            sp <- split(v, control$strata)
            res <- matrix(nrow = nperms, ncol = n)
            for(i in seq_len(perm.wi))
                res[i,] <- sapply(sp, function(x) x[ord[i,]])
        } else {
            ## different permutations within blocks
            tab <- table(control$strata)
            ng <- length(tab)
            pg <- unique(tab)
            if(length(pg) > 1) {
                ## different number of observations per level of strata
                if(control$type == "grid")
                    ## FIXME: this should not be needed once all checks are
                    ## in place in permCheck()
                    stop("Unbalanced grid designs are not supported")
                control.wi <- permControl(type = control$type,
                                          mirror = control$mirror)
                sp <- split(v, control$strata)
                res <- vector(mode = "list", length = ng)
                add <- c(0, cumsum(tab)[1:(ng-1)])
                for(j in seq(along = tab)) {
                    ord <- switch(control.wi$type,
                                  free = all.free(tab[j]),
                                  series = all.series(tab[j],
                                  control=control.wi))
                    perm.wi <- nrow(ord)
                    if(j == 1) {
                        a <- 1
                        b <- nperms / perm.wi
                    } else {
                        b <- b/perm.wi
                        a <- nperms / (b*perm.wi)
                    }
                    res[[j]] <- matrix(rep(bar(ord+add[j], a),
                                           each = b),
                                       ncol = tab[j])
                }
                res <- do.call(cbind, res)
            } else {
                ## same number of observations per level of strata
                control.wi <- permControl(type = control$type,
                                          mirror = control$mirror,
                                          nrow = control$nrow,
                                          ncol = control$ncol)
                ord <- switch(control$type,
                              free = all.free(pg),
                              series = all.series(pg, control = control.wi),
                              grid = all.grid(pg, control = control.wi)
                              )
                perm.wi <- nrow(ord)
                add <- seq(from = 0, by = pg, length.out = ng)
                res <- vector(mode = "list", length = ng)
                a <- 1
                b <- nperms / perm.wi
                for(i in seq_len(ng)) {
                    res[[i]] <- matrix(rep(bar(ord+add[i], a), each = b),
                                       ncol = pg)
                    a <- a*perm.wi
                    b <- b/perm.wi
                }
                res <- do.call(cbind, res)
            }
        }
    } else {
        ## not permuting within blocks or are permuting strata
        res <- switch(type,
                      free = all.free(n),
                      series = all.series(n, control=control),
                      grid = all.grid(n, control=control),
                      strata = all.strata(n, control=control)
                      )
    }
    ## some times storage.mode of res is numeric, sometimes
    ## it is integer, set to "integer" for comparisons using
    ## identical to match the observed ordering
    storage.mode(res) <- "integer"
    if(!observed) {
        obs.row <- apply(res, 1, function(x, v) {identical(x, v)}, v)
        res <- res[!obs.row, ]
        ## reduce the number of permutations to get rid of the
        ## observed ordering
        control$nperm <- control$nperm - 1
    }
    class(res) <- "allPerms"
    attr(res, "control") <- control
    attr(res, "observed") <- observed
    res
}

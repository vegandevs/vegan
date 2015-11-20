`mrpp` <-
    function (dat, grouping, permutations = 999, distance = "euclidean", 
              weight.type = 1, strata = NULL,
              parallel = getOption("mc.cores")) 
{
    EPS <- sqrt(.Machine$double.eps)
    classmean <- function(ind, dmat, indls) {
        sapply(indls, function(x)
               mean(c(dmat[ind == x, ind == x]),
                    na.rm = TRUE))
    }
    mrpp.perms <- function(ind, dmat, indls, w) {
        weighted.mean(classmean(ind, dmat, indls), w = w, na.rm = TRUE)
    }
    if (inherits(dat, "dist")) 
        dmat <- dat
    else if ((is.matrix(dat) || is.data.frame(dat)) &&
               isSymmetric(unname(as.matrix(dat)))) {
        dmat <- dat
        attr(dmat, "method") <- "user supplied square matrix"
    }
    else dmat <- vegdist(dat, method = distance)
    if (any(dmat < -sqrt(.Machine$double.eps)))
        stop("dissimilarities must be non-negative")
    distance <- attr(dmat, "method")
    dmat <- as.matrix(dmat)
    diag(dmat) <- NA
    N <- nrow(dmat)
    grouping <- factor(grouping)
    indls <- levels(grouping)
    ncl <- sapply(indls, function(x) sum(grouping == x))
    w <- switch(weight.type, ncl, ncl - 1, ncl * (ncl - 1)/2)
    classdel <- classmean(grouping, dmat, indls)
    names(classdel) <- names(ncl) <- indls
    del <- weighted.mean(classdel, w = w, na.rm = TRUE)
    E.del <- mean(dmat, na.rm = TRUE)
    ## 'Classification strength' if weight.type == 1
    ## Do not calculate classification strength because there is no
    ## significance test for it. Keep the item in reserve for
    ## possible later re-inclusion.
    CS <- NA
    permutations <- getPermuteMatrix(permutations, N, strata = strata)
    if (ncol(permutations) != N)
        stop(gettextf("'permutations' have %d columns, but data have %d rows",
                      ncol(permutations), N))

    control <- attr(permutations, "control")
    if(nrow(permutations)) {
        perms <- apply(permutations, 1, function(indx) grouping[indx])
        permutations <- ncol(perms)

        ## Parallel processing
        if (is.null(parallel))
            parallel <- 1
        hasClus <- inherits(parallel, "cluster")
        if (hasClus || parallel > 1) {
            if(.Platform$OS.type == "unix" && !hasClus) {
                m.ds <- unlist(mclapply(1:permutations, function(i, ...)
                                        mrpp.perms(perms[,i], dmat, indls, w),
                                        mc.cores = parallel))
            } else {
                if (!hasClus) {
                    parallel <- makeCluster(parallel)
                }
                m.ds <- parCapply(parallel, perms, function(x)
                                  mrpp.perms(x, dmat, indls, w))
                if (!hasClus)
                    stopCluster(parallel)
            }
        } else {
            m.ds <- apply(perms, 2, function(x) mrpp.perms(x, dmat, indls, w))
        }
        p <- (1 + sum(del + EPS >= m.ds))/(permutations + 1)
        r2 <- 1 - del/E.del
    } else { # no permutations
        m.ds <- p <- r2 <- NA
        permutations <- 0
    }
    out <- list(call = match.call(), delta = del, E.delta = E.del, CS = CS,
                n = ncl, classdelta = classdel, Pvalue = p, A = r2,
                distance = distance, weight.type = weight.type,
                boot.deltas = m.ds, permutations = permutations,
                control = control)
    class(out) <- "mrpp"
    out
}

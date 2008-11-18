"mrpp" <-
function (dat, grouping, permutations = 1000, distance = "euclidean", 
    weight.type = 1, strata) 
{
    mrpp.perms <- function(ind, dmat, indls, w) {
        weighted.mean(sapply(indls, function(x) mean(c(dmat[ind == 
            x, ind == x]), na.rm = TRUE)), w = w, na.rm = TRUE)
    }
    if (inherits(dat, "dist")) 
        dmat <- dat
    else if (is.matrix(dat) && nrow(dat) == ncol(dat) && all(dat[lower.tri(dat)] == 
        t(dat)[lower.tri(dat)])) {
        dmat <- dat
        attr(dmat, "method") <- "user supplied square matrix"
    }
    else dmat <- vegdist(dat, method = distance)
    distance <- attr(dmat, "method")
    dmat <- as.matrix(dmat)
    diag(dmat) <- NA
    N <- nrow(dmat)
    ind <- as.numeric(grouping)
    indls <- unique(ind)
    w <- sapply(indls, function(x) sum(ind == x))
    w <- switch(weight.type, w, w - 1, w * (w - 1)/2)
    del <- mrpp.perms(ind, dmat, indls, w)
    if (missing(strata)) 
        strata <- NULL
    perms <- sapply(1:permutations, function(x) ind[permuted.index(N, 
        strata = strata)])
    m.ds <- numeric(permutations)
    m.ds <- apply(perms, 2, function(x) mrpp.perms(x, dmat, indls, 
        w))
    E.del <- mean(m.ds)
    p <- (1 + sum(del >= m.ds))/(permutations + 1)
    r2 <- 1 - del/E.del
    out <- list(call = match.call(), delta = del, E.delta = E.del, 
        Pvalue = p, A = r2, distance = distance, weight.type = weight.type, 
        boot.deltas = m.ds, permutations = permutations)
    if (!is.null(strata)) {
        out$strata <- deparse(substitute(strata))
        out$stratum.values <- strata
    }
    class(out) <- "mrpp"
    out
}


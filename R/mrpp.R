"mrpp" <-
function (dat, grouping, permutations = 1000, distance = "euclidean", 
    weight.type = 1, strata) 
{
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
    grouping <- factor(grouping)
    indls <- levels(grouping)
    ncl <- sapply(indls, function(x) sum(grouping == x))
    w <- switch(weight.type, ncl, ncl - 1, ncl * (ncl - 1)/2)
    classdel <- classmean(grouping, dmat, indls)
    names(classdel) <- names(ncl) <- indls
    del <- weighted.mean(classdel, w = w, na.rm = TRUE)
    E.del <- mean(dmat, na.rm = TRUE)
    ## 'Classification strength' if weight.type == 3
    if (weight.type == 3) {
        CS <- N*(N-1)/2*(E.del - del)/(N*(N-1)/2 - sum(w))
    } else {
        CS <- NA
    }
    if (missing(strata)) 
        strata <- NULL
    perms <- sapply(1:permutations, function(x) grouping[permuted.index(N, 
        strata = strata)])
    m.ds <- numeric(permutations)
    m.ds <- apply(perms, 2, function(x) mrpp.perms(x, dmat, indls, 
        w))
    p <- (1 + sum(del >= m.ds))/(permutations + 1)
    r2 <- 1 - del/E.del
    out <- list(call = match.call(), delta = del, E.delta = E.del, CS = CS,
        n = ncl, classdelta = classdel,
        Pvalue = p, A = r2, distance = distance, weight.type = weight.type, 
        boot.deltas = m.ds, permutations = permutations)
    if (!is.null(strata)) {
        out$strata <- deparse(substitute(strata))
        out$stratum.values <- strata
    }
    class(out) <- "mrpp"
    out
}

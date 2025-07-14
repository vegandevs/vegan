`vectorfit` <-
    function (X, P, permutations = 0, strata = NULL, w, ...)
{
    EPS <- sqrt(.Machine$double.eps)
    if (missing(w) || is.null(w))
        w <- 1
    if (length(w) == 1)
        w <- rep(w, nrow(X))
    P <- as.matrix(P)
    if (nrow(P) != nrow(X))
        stop("input data have non-matching numbers of observations")
    Xw <- .Call(do_wcentre, X, w, PACKAGE = "vegan")
    Pw <- .Call(do_wcentre, P, w, PACKAGE = "vegan")
    colnames(Pw) <- colnames(P)
    nc <- ncol(X)
    Q <- qr(Xw)
    H <- qr.fitted(Q, Pw)
    heads <- qr.coef(Q, Pw)
    r <- diag(cor(H, Pw)^2)
    r[is.na(r)] <- 0
    heads <- decostand(heads, "norm", 2)
    heads <- t(heads)
    if (is.null(colnames(X)))
        colnames(heads) <- paste("Dim", 1:nc, sep = "")
    else colnames(heads) <- colnames(X)
    ## make permutation matrix for all variables handled in the next loop
    nr <- nrow(X)
    permat <- getPermuteMatrix(permutations, nr, strata = strata)
    if (ncol(permat) != nr)
        stop(gettextf("'permutations' have %d columns, but data have %d rows",
                          ncol(permat), nr))
    permutations <- nrow(permat)

    if (permutations) {
        ptest <- function(indx, ...) {
            take <- P[indx, , drop = FALSE]
            take <- .Call(do_wcentre, take, w, PACKAGE = "vegan")
            Hperm <- qr.fitted(Q, take)
            diag(cor(Hperm, take))^2
        }
        permstore <- sapply(1:permutations, function(indx, ...) ptest(permat[indx,], ...))
        ## Single variable is dropped to a vector, and otherwise
        ## permutations are the matrix columns and variables are rows
        if (!is.matrix(permstore))
            permstore <- matrix(permstore, ncol=permutations)
        permstore <- sweep(permstore, 1, r - EPS, ">=")
        validn <- rowSums(is.finite(permstore))
        pvals <- (rowSums(permstore, na.rm = TRUE) + 1)/(validn + 1)
    }
    else pvals <- NULL
    sol <- list(arrows = heads, r = r, permutations = permutations,
                pvals = pvals)
    sol$control <- attr(permat, "control")
    class(sol) <- "vectorfit"
    sol
}

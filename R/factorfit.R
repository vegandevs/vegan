`factorfit` <-
    function (X, P, permutations = 0, strata = NULL, w,  ...)
{
    EPS <- sqrt(.Machine$double.eps)
    P <- as.data.frame(P)
    ## Check that all variables are factors, and coerce if necessary
    if(any(!sapply(P, is.factor)))
        P <- data.frame(lapply(P, function(x)
                        if (is.factor(x)) x else factor(x)))
    P <- droplevels(P, exclude = NA) ## make sure only the used levels are present
    if (any(!sapply(P, is.factor)))
        stop("all non-numeric variables must be factors")
    NR <- nrow(X)
    NC <- ncol(X)
    NF <- ncol(P)
    if (missing(w) || is.null(w))
        w <- 1
    if (length(w) == 1)
        w <- rep(w, NR)
    r <- NULL
    pval <- NULL
    totvar <- .Call(do_goffactor, X, rep(1L, NR), 1L, w)
    sol <- centroids.cca(X, P, w)
    var.id <- rep(names(P), sapply(P, nlevels))
    ## make permutation matrix for all variables handled in the next loop
    permat <- getPermuteMatrix(permutations, NR, strata = strata)
    permutations <- nrow(permat)

    for (i in seq_along(P)) {
        A <- as.integer(P[[i]])
        NL <- nlevels(P[[i]])
        invar <- .Call(do_goffactor, X, A, NL, w)
        r.this <- 1 - invar/totvar
        r <- c(r, r.this)
        if (permutations) {
            A <- as.integer(P[[i]])
            NL <- nlevels(P[[i]])
            ptest <- function(indx, ...) {
                take <- A[indx]
                invar <- .Call(do_goffactor, X, take, NL, w)
                1 - invar/totvar
            }
            tmp <- sapply(seq_len(permutations),
                          function(indx,...) ptest(permat[indx,], ...))
            pval.this <- (sum(tmp >= r.this - EPS) + 1)/(permutations + 1)
            pval <- c(pval, pval.this)
        }
    }
    if (is.null(colnames(X)))
        colnames(sol) <- paste("Dim", 1:ncol(sol), sep = "")
    else colnames(sol) <- colnames(X)
    names(r) <- names(P)
    if (!is.null(pval))
        names(pval) <- names(P)
    out <- list(centroids = sol, r = r, permutations = permutations,
                pvals = pval, var.id = var.id)
    out$control <- attr(permat, "control")
    class(out) <- "factorfit"
    out
}

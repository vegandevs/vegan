"factorfit" <-
    function (X, P, permutations = 0, strata, choices = c(1, 2),
              display = c("sites","lc"), w = weights(X),  ...) 
{
    weights.default <- function(object, ...) NULL
    display <- match.arg(display)
    w <- eval(w)
    P <- as.data.frame(P)
    if (any(!sapply(P, is.factor))) 
        stop("All fitted variables must be factors")
    X <- scores(X, display = display, choices, ...)
    NR <- nrow(X)
    NC <- ncol(X)
    NF <- ncol(P)
    if (is.null(w))
        w <- 1
    if (length(w) == 1)
        w <- rep(w, NR)
    r <- NULL
    pval <- NULL
    totvar <- .C("goffactor", as.double(X), as.integer(rep(0, 
                                                           NR)), as.double(w), as.integer(NR), as.integer(NC), as.integer(1), 
                 double(1), double(1), double(1), var = double(1), PACKAGE = "vegan")$var
    sol <- centroids.cca(X, P, w)
    var.id <- rep(names(P), sapply(P, nlevels))
    for (i in 1:length(P)) {
        A <- as.integer(P[[i]])
        NL <- nlevels(P[[i]])
        invar <- .C("goffactor", as.double(X), as.integer(A - 
                                                          1), as.double(w), as.integer(NR), as.integer(NC), 
                    as.integer(NL), double(NL), double(NL), double(NL), 
                    var = double(1), PACKAGE = "vegan")$var
        r.this <- 1 - invar/totvar
        r <- c(r, r.this)
        if (permutations) {
            A <- as.integer(P[[i]])
            NL <- nlevels(P[[i]])
            tmp <- rep(NA, permutations)
            for (i in 1:permutations) {
                indx <- permuted.index(length(A), strata)
                take <- A[indx]
                invar <- .C("goffactor", as.double(X), as.integer(take - 
                                                                  1), as.double(w), as.integer(NR), as.integer(NC), 
                            as.integer(NL), double(NL), double(NL), double(NL), 
                            var = double(1), PACKAGE = "vegan")$var
                tmp[i] <- 1 - invar/totvar
            }
            pval.this <- (sum(tmp > r.this) + 1)/(permutations + 1)
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
    if (!missing(strata)) {
        out$strata <- deparse(substitute(strata))
        out$stratum.values <- strata
    }
    class(out) <- "factorfit"
    out
}

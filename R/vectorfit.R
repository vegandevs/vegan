"vectorfit" <-
    function (X, P, permutations = 0, strata, choices = c(1, 2), 
              display = c("sites", "lc"), w = weights(X), ...) 
{
    weights.default <- function(object, ...) NULL
    display <- match.arg(display)
    w <- eval(w)
    X <- scores(X, display = display, choices, ...)
    if (is.null(w)) 
        w <- 1
    if (length(w) == 1) 
        w <- rep(w, nrow(X))
    Xw <- .C("wcentre", x = as.double(X), as.double(w), as.integer(nrow(X)),
             as.integer(ncol(X)), PACKAGE = "vegan")$x
    dim(Xw) <- dim(X)
    P <- as.matrix(P)
    Pw <- .C("wcentre", x = as.double(P), as.double(w), as.integer(nrow(P)),
             as.integer(ncol(P)), PACKAGE = "vegan")$x
    dim(Pw) <- dim(P)
    colnames(Pw) <- colnames(P)
    nc <- ncol(X)
    Q <- qr(Xw)
    H <- qr.fitted(Q, Pw)
    heads <- qr.coef(Q, Pw)
    r <- diag(cor(H, Pw)^2)
    heads <- decostand(heads, "norm", 2)
    heads <- t(heads)
    if (is.null(colnames(X))) 
        colnames(heads) <- paste("Dim", 1:nc, sep = "")
    else colnames(heads) <- colnames(X)
    if (permutations) {
        nr <- nrow(X)
        permstore <- matrix(nrow = permutations, ncol = ncol(P))
        for (i in 1:permutations) {
            indx <- permuted.index(nrow(P), strata)
            take <- P[indx, , drop = FALSE]
            take <- .C("wcentre", x = as.double(take), as.double(w),
                       as.integer(nrow(take)), as.integer(ncol(take)),
                       PACKAGE = "vegan")$x
            dim(take) <- dim(P)
            Hperm <- qr.fitted(Q, take)
            permstore[i, ] <- diag(cor(Hperm, take))^2
        }
        permstore <- sweep(permstore, 2, r, ">")
        pvals <- (apply(permstore, 2, sum) + 1)/(permutations + 1)
    }
    else pvals <- NULL
    sol <- list(arrows = heads, r = r, permutations = permutations, 
                pvals = pvals)
    if (!missing(strata)) {
        sol$strata <- deparse(substitute(strata))
        sol$stratum.values <- strata
    }
    class(sol) <- "vectorfit"
    sol
}

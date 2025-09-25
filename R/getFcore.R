### getF.c: prototype of essential commands in R

#' R prototypes of (alternative) implementations of getF.c
#'
#' Function `getFcore` is the R prototype of getF.c as implemented in
#' vegan release 2.6-6. This version permutes together response Y,
#' conditions Z and weights w, but keeps constraints X
#' non-permuted. However, X must be reweighted by shuffled w and
#' therefore we need a new QR decomposition of X. This is only needed
#' in weighted ordination (CCA) and in RDA we can re-use the same QR.
#'
#' Function `XgetFcore` implements an alternative scheme, where Y, Z &
#' w are non-shuffled, and only X is shuffled. The shuffled X must be
#' re-weighted by non-shuffled weights w. Because X is shuffled, we
#' need a new QR decomposition also in unweighted analysis (RDA,
#' dbRDA).
#'
#' With same permutations, functions `anova.cca`, `permutest.cca`,
#' `getFcore` and `getFcore` return same permutation F values.
#'
#' @examples
#' library(vegan)
#' data(mite, mite.env)
#' perm <- shuffleSet(nrow(mite), 999)
#' mod <- cca(mite ~ SubsDens + WatrCont + Condition(Topo + Shrub),
#'    data=mite.env)
#' ano0 <- anova(mod, permutations=perm)
#' ano0 <- drop(permustats(mod)$permutations)
#' anoYZw <- getFcore(mod, perm)
#' anoX <- XgetFcore(mod, perm)
#' plot(ano0, anoYXw); abline(0,1)
#' plot(ano0, anoX); abline(0,1)
#' ## evenness of row weights
#' diversity(rowSums(mite), "inv") # virtual N
#' nrow(mite) # N in unweighted model
#'
#' @param m fitted partial constrained ordination model
#' @param p permutation matrix

`getFcore` <-
    function(m, p)
{
    ## data
    Y <- ordiYbar(m, "partial") # reduced model
    QZ <- m$pCCA$QR
    QR <- m$CCA$QR
    w <- weights(m)

    ## permutations
    if (missing(p))
        p <- matrix(seq_len(nrow(Y)), nrow = 1)
    niter <- nrow(p)
    ss <- numeric(niter)

    ## Set up before the loop
    Z <- qr.X(QZ) # weighted Z
    X <- .Call(test_qrXw, QR, w, ncol(Z)) # unweighted [X]

    for (iter in seq_len(niter)) {
        ## permute Y & w
        Yperm <- Y[p[iter,],]
        Zperm <- Z[p[iter,],, drop = FALSE]
        wperm <- w[p[iter,]]
        ## Partial
        QZ <- qr(Zperm)
        Yperm <- qr.resid(QZ, Yperm)
        ## Constrained
        Xrew <- .Call(do_wcentre, X, wperm)
        Xrew <- qr.resid(QZ, Xrew) # "residualized predictor" X
        QR <- qr(Xrew)
        Yfit <- qr.fitted(QR, Yperm)
        Yres <- qr.resid(QR, Yperm)
        ss[iter] <- sum(Yfit^2)/sum(Yres^2)
    }
    df <- (nobs(m) - m$CCA$QR$rank - 1) / m$CCA$rank
    ss * df
}

`XgetFcore` <-
    function(m, p)
{
    ## data
    Y <- ordiYbar(m, "partial") # reduced model
    QZ <- m$pCCA$QR
    QR <- m$CCA$QR
    w <- weights(m)
    ## permutations
    if (missing(p))
        p <- matrix(seq_len(nrow(Y)), nrow = 1)
    niter <- nrow(p)
    ss <- numeric(niter)

    ## Setting up before the loop
    X <- .Call(test_qrXw, QR, w, m$pCCA$rank)

    for (iter in seq_len(niter)) {
        ## permute X: order() makes X-permutation equal to
        ## Y-permutation, and order(order(i)) == i
        Xperm <- X[order(p[iter,]),]
        Xrew <- .Call(do_wcentre, Xperm, w)
        Xrew <- qr.resid(QZ, Xrew)
        QR <- qr(Xrew)
        Yfit <- qr.fitted(QR, Y)
        Yres <- qr.resid(QR, Y)
        ss[iter] <- sum(Yfit^2)/sum(Yres^2)
    }
    df <- (nobs(m) - m$CCA$QR$rank - 1) / m$CCA$rank
    ss * df
}

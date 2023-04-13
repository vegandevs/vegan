## getF.c: essential commands in R

#' Inspect changing getF.c
#'
#' This version permutes Y & w, does not permute Z, but reweights Z &
#' X. This is the first version in this branch that fixes problems
#' with weights in most cases, but fails in the most extreme of Cajo's
#' tests (Pinho). This is point-to-point identical to C code as
#' implemement in branch biased-anova-cca v2.6-3-29-gaaf6f700. This is
#' similar to vegan tests prior to 2.5-1.
#'
#' @examples
#' library(vegan)
#' data(mite, mite.env)
#' perm <- shuffleSet(nrow(mite), 999)
#' mod <- cca(mite ~ SubsDens + WatrCont + Condition(Topo + Shrub),
#'    data=mite.env)
#' getFcore(mod, perm)
#' ## evenness of row weights
#' diversity(rowSums(mite), "inv") # virtual N
#' nrow(mite) # N in unweighted model
#'
#' @param m fitted constrained ordination model
#' @param p permutation matrix

`getFcore` <-
    function(m, p)
{
    ## data
    Y <- ordiYbar(m, "partial") # reduced model
    QZ <- m$pCCA$QR
    QR <- m$CCA$QR
    w <- weights(m)

    ## Set up before the loop
    ## Z <- .Call(vegan:::test_qrXw, QZ, w, 0)
    ## X <- .Call(vegan:::test_qrXw, QR, w, ncol(Z))

    ## permutations
    if (missing(p))
        p <- matrix(seq_len(nrow(Y)), nrow = 1)
    niter <- nrow(p)
    ss <- numeric(niter)

    ## Set up before the loop
    Z <- qr.X(QZ) # weighted Z
    X <- .Call(test_qrXw, QR, w, ncol(Z)) # unweighted [ZX]

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
    list(P = (sum(ss >= m$CCA$tot.chi/m$CA$tot.chi) + 1) / (niter + 1),
         ss = ss)
}

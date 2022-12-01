## getF.c: essential commands in R

#' Inspect changing getF.c
#'
#' This version permutes Y but does not permute Z or reweight Z &
#' X. This is the original implementation of getF.c in vegan since
#' 2.5-1.
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
    ## permutations
    if (missing(p))
        p <- matrix(seq_len(nrow(Y)), nrow = 1)
    niter <- nrow(p)
    ss <- numeric(niter)

    for (iter in seq_len(niter)) {
        ## permute Y, Z & w
        Yperm <- Y[p[iter,],]
        Yperm <- qr.resid(QZ, Yperm)
        Yfit <- qr.fitted(QR, Yperm)
        Yres <- qr.resid(QR, Yperm)
        ss[iter] <- sum(Yfit^2)/sum(Yres^2)
    }
    list(P = (sum(ss >= m$CCA$tot.chi/m$CA$tot.chi) + 1) / (niter + 1),
         ss = ss)
}

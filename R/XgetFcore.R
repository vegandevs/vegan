## getF.c: essential commands in R

#' Inspect changing getF.c
#'
#' This version permutes and reweights constraints X.
#'
#' @param m fitted constrained ordination model
#' @param p permutation matrix

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
    ## Z <- qr.X(QZ) ##.Call(vegan:::test_qrXw, QZ, w) # not needed
    X <- .Call(vegan:::test_qrXw, QR, w, QZ$rank)

    for (iter in seq_len(niter)) {
        ## permute X: order() makes X-permutation equal to
        ## Y-permutation, and order(order(i)) == i
        Xperm <- X[order(p[iter,]),]
        ## The following is done in a loop which reorders rows of X
        ##Zrew <- .Call(vegan:::do_wcentre, Z, wperm)
        ## QZ <- qr(Zperm) # Z will be constant when permuting X
        # Yperm <- qr.resid(QZ, Yperm) #
        Xrew <- .Call(vegan:::do_wcentre, Xperm, w)
        Xrew <- qr.resid(QZ, Xrew)
        QR <- qr(Xrew)
        Yfit <- qr.fitted(QR, Y)
        Yres <- qr.resid(QR, Y)
        ss[iter] <- sum(Yfit^2)/sum(Yres^2)
    }
    list(P = (sum(ss >= m$CCA$tot.chi/m$CA$tot.chi) + 1) / (niter + 1),
         ss = ss)
}

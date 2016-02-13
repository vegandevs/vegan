`ade2vegancca` <-
    function(object)
{
    nf <- object$nf
    CCA <- list(eig = object$eig,
                u = as.matrix(object$l1),
                v = as.matrix(object$c1),
                wa = sweep(as.matrix(object$ls), 2,
                           1/sqrt(object$eig[1:nf]), "*"),
                biplot = as.matrix(object$cor)[-1,],
                rank = object$rank,
                tot.chi = sum(object$eig),
                QR = NA,
                envcentre = NA,
                Xbar = NA)
    out <- list(call = object$call,
                grand.total = NA,
                rowsum = object$lw,
                colsum = object$cw,
                tot.chi = NA,
                pCCA = NULL,
                CCA = CCA,
                CA = NULL,
                method = "cca",
                inertia = "mean squared contingency coefficient")
    class(out) = c("ade4cca", "cca")
    out
}

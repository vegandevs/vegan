`summary.cca` <-
    function (object, scaling = 2, axes = 6, display=c("sp","wa","lc","bp","cn"), 
              digits = max(3, getOption("digits") - 3), ...) 
{
    if (inherits(object, "pcaiv")) {
        warning("this is an ade4 object which vegan cannot handle")
        axes <- min(axes, object$nf)
        object <- ade2vegancca(object)
    }
    axes <- min(axes, sum(object$CCA$rank, object$CA$rank))
    summ <- list()
    if (axes && length(display) && (!is.na(display) && !is.null(display))) 
        summ <- scores(object, scaling = scaling, choices = 1:axes, display = display, ...)
    ## scores() drops list to a matrix if there is only one item: workaround below.
    if (!is.list(summ) && length(display) == 1) {
        nms <- c("species", "sites", "constraints", "biplot", "centroids")
        names(nms) <- c("sp","wa","lc","bp","cn")
        summ <- list(summ)
        names(summ) <- nms[display]
    }
    summ$call <- object$call
    summ$tot.chi <- object$tot.chi
    summ$partial.chi <- object$pCCA$tot.chi
    summ$constr.chi <- object$CCA$tot.chi
    summ$unconst.chi <- object$CA$tot.chi
    summ$ev.con <- object$CCA$eig
    summ$ev.uncon <- object$CA$eig
    ev.account <- summ$tot.chi
    if (!is.null(object$pCCA)) 
        ev.account <- ev.account - summ$partial.chi
    if (!is.null(object$CCA))
        summ$ev.con.account <- cumsum(summ$ev.con)/ev.account
    summ$ev.uncon.account <-
        (max(summ$constr.chi, 0) + cumsum(summ$ev.uncon))/ev.account
    if (!is.null(object$CCA))
        summ$cca.acc <- cumsum(summ$ev.con)/summ$constr.chi
    summ$ev.head <- c(summ$ev.con, summ$ev.uncon)[1:axes]
    summ$scaling <- scaling
    summ$digits <- digits
    summ$inertia <- object$inertia
    summ$method <- object$method
    class(summ) <- "summary.cca"
    summ
}

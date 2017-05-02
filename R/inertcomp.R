`inertcomp` <-
    function (object, display = c("species", "sites"),
              statistic = c("explained", "distance"), proportional = FALSE)
{
    display <- match.arg(display)
    statistic <- match.arg(statistic)
    if (!inherits(object, "cca"))
        stop("can be used only with objects inheriting from 'cca'")
    if (inherits(object, c("capscale", "dbrda")) && display == "species")
        stop(gettextf("cannot analyse species with '%s'", object$method))
    what <- if(display == "species") "v" else "u"
    w <- weights(object, display = display)
    pCCA <- object$pCCA$Fit
    CCA <- object$CCA[[what]]
    CA <- object$CA[[what]]
    ## row names will be lost later: save here
    labels <- if(!is.null(CA))
                  rownames(CCA)
              else if(!is.null(CCA))
                  rownames(CA)
              else
                  rownames(pCCA)
    ## imaginary dimensions for dbrda
    if (inherits(object, "dbrda")) {
        CCA <- cbind(CCA, object$CCA$imaginary.u)
        CA <- cbind(CA, object$CA$imaginary.u)
    }
    if (inherits(object, "rda")) {
        nr <- nobs(object) - 1
    }
    else {
        nr <- 1
    }
    if (!is.null(pCCA) && object$pCCA$rank > 0) {
        if (display == "sites")
            pCCA <- t(pCCA)
        if (inherits(object, "dbrda"))
            pCCA <- diag(pCCA)
        else
            pCCA <- diag(crossprod(pCCA))/nr
    }
    if (!is.null(CCA) && object$CCA$rank > 0)
        CCA <- rowSums(diag(w, length(w)) %*% CCA^2 %*%
                           diag(object$CCA$eig, length(object$CCA$eig)))
    if (!is.null(CA) && object$CA$rank > 0)
        CA <- rowSums(diag(w, length(w)) %*% CA^2 %*%
                          diag(object$CA$eig, length(object$CA$eig)))
    out <- cbind(pCCA, CCA, CA)
    if (statistic == "distance" && !proportional) {
        w <- weights(object, display = display)
        if (display == "sites" &&
            !is.null(object$na.action) &&
            inherits(object$na.action, "exclude"))
            w <- w[-object$na.action]
        out <- sweep(out, 1, w, "/")
    }
    if (proportional)
        out <- sweep(out, 1, rowSums(out), "/")
    ## get back names
    rownames(out) <- labels
    out
}

"inertcomp" <-
    function (object, display = c("species", "sites"), statistic = c("explained", 
                                                       "distance"), proportional = FALSE) 
{
    display <- match.arg(display)
    statistic <- match.arg(statistic)
    if (!inherits(object, "cca"))
        stop("can be used only with objects inheriting from 'cca'")
    if (inherits(object, "capscale") && display == "species")
        stop("cannot analyse species with 'capscale'")
    pCCA <- object$pCCA$Fit
    CCA <- object$CCA$Xbar
    CA <- object$CA$Xbar
    if (inherits(object, "rda")) {
        nr <- nrow(CA) - 1
        if (is.null(nr)) 
            nr <- nrow(CCA) - 1
        if (is.null(nr)) 
            nr <- nrow(pCCA) - 1
    }
    else {
        nr <- 1
    }
    if (!is.null(pCCA)) {
        if (display == "sites") 
            pCCA <- t(pCCA)
        pCCA <- diag(crossprod(pCCA))/nr
    }
    if (!is.null(CCA)) {
        CCA <- qr.fitted(object$CCA$QR, CCA)
        if (display == "sites") 
            CCA <- t(CCA)
        CCA <- diag(crossprod(CCA))/nr
    }
    if (!is.null(CA)) {
        if (display == "sites") 
            CA <- t(CA)
        CA <- diag(crossprod(CA))/nr
    }
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
    out
}

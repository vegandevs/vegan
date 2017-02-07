`inertcomp` <-
    function (object, display = c("species", "sites"),
              unity = FALSE, proportional = FALSE)
{
    display <- match.arg(display)
    if (!inherits(object, "cca"))
        stop("can be used only with objects inheriting from 'cca'")
    if (inherits(object, c("capscale", "dbrda")) && display == "species")
        stop(gettextf("cannot analyse species with '%s'", object$method))
    if (inherits(object, "dbrda"))
        display <- "dbrda"
    ## function to get the eigenvalues
    getComps <- function(x, display) {
        if(!is.null(x))
            switch(display,
                   "species" = colSums(x^2),
                   "sites" = rowSums(x^2),
                   "dbrda" = diag(x)
                   )
    }
    pCCA <- object$pCCA$Fit
    tot <- if (!is.null(object$CCA)) object$CCA$Xbar else object$CA$Xbar
    if (!is.null(pCCA))
        tot <- tot + pCCA
    CCA <- qr.fitted(object$CCA$QR, object$CCA$Xbar)
    if (inherits(object, "dbrda"))
        CCA <- qr.fitted(object$CCA$QR, t(CCA))
    CA <- object$CA$Xbar

    out <- cbind("pCCA" = getComps(pCCA, display),
                 "CCA" = getComps(CCA, display),
                 "CA" = getComps(CA, display))
    if (unity)
        out <- out/sum(out)
    if (proportional)
        out <- out/getComps(tot, display)
    out
}

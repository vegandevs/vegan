`goodness.rda` <-
    function (object, display = c("species", "sites"), choices, model = c("CCA", 
                                                                "CA"), statistic = c("explained", "distance"), summarize = FALSE, 
              ...) 
{
    model <- match.arg(model)
    display <- match.arg(display)
    if (inherits(object, "capscale") && display == "species") 
        stop("display = \"species\" not available for 'capscale'")
    if (is.null(object$CCA)) 
        model <- "CA"
    if (is.null(object[[model]])) 
        stop("model ", model, " is not available")
    statistic <- match.arg(statistic)
    cs <- weights(object, display = display)
    if (display == "species") {
        if (is.null(object$CCA)) 
            Xbar <- object$CA$Xbar
        else Xbar <- object$CCA$Xbar
        v <- object[[model]]$v.eig
        tot <- diag(crossprod(Xbar)/(nrow(Xbar) - 1))
    }
    else {
        Xbar <- object$CA$Xbar
        tot <- diag(crossprod(t(Xbar)))
        if (!is.null(tot)) 
            tot <- tot/(nrow(Xbar) - 1)
        if (!is.null(object$CCA)) {
            Xbar <- object$CCA$Xbar
            Xbar <- qr.fitted(object$CCA$QR, Xbar)
            tot <- tot + diag(crossprod(t(Xbar)))/(nrow(Xbar) - 
                                                   1)
        }
        v <- object[[model]]$u.eig
    }
    if (!missing(choices)) 
        v <- v[, choices, drop = FALSE]
    vexp <- t(apply(v^2, 1, cumsum))
    vexp <- sweep(vexp, 1, cs, "*")
    if (statistic == "explained") {
        if (!is.null(object$pCCA)) {
            Xbar <- object$pCCA$Fit
            if (display == "sites") 
                Xbar <- t(Xbar)
            ptot <- diag(crossprod(Xbar))
            tot <- tot + ptot
            vexp <- sweep(vexp, 1, ptot, "+")
        }
        vexp <- sweep(vexp, 1, tot, "/")
    }
    else {
        if (display == "sites" && (!is.null(object$CCA) || !is.null(object$pCCA)))
            stop("statistic 'distance' not available for sites in constrained analysis")
        vexp <- sweep(-(vexp), 1, tot, "+")
        vexp[vexp < 0] <- 0
        vexp <- sweep(sqrt(vexp), 1, cs, "/")
    }
    if (summarize) 
        vexp <- vexp[, ncol(vexp)]
    vexp
}


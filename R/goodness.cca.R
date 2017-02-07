`goodness.cca` <-
    function (object, choices, display = c("species", "sites"),
              model = c("CCA", "CA"),
              summarize = FALSE, addprevious = FALSE, ...)
{
    display <- match.arg(display)
    model <- match.arg(model)
    if (!inherits(object, "cca"))
        stop("can be used only with objects inheriting from 'cca'")
    if (inherits(object, c("capscale", "dbrda")) && display == "species")
        stop(gettextf("cannot analyse species with '%s'", object$method))
    v <- sqrt(weights(object, display="species")) * object[[model]]$v
    if (is.null(v))
        stop(gettextf("model = '%s' deos not exist", model))
    if (display == "sites")
        u <- sqrt(weights(object, display="sites")) * object[[model]]$u
    eig <- object[[model]]$eig
    if (!inherits(object, "dbrda"))
        eig <- eig[eig > 0]
    ## imaginary dimensions for dbrda
    if (inherits(object, "dbrda"))
        CA <- cbind(CA, object[[model]][["imaginary.u"]])
    ## take only chosen axes within the component
    if (!missing(choices)) {
        choices <- choices[choices <= ncol(CA)]
        v <- v[, choices, drop = FALSE]
        if (display == "sites")
            u <- u[, choices, drop = FALSE]
        eig <- eig[choices]
    }
    if (addprevious) {
        if (!is.null(object$pCCA))
            prev <- prev$pCCA$Fit
        else
            prev <- 0
        if (model == "CA" && !is.null(object$CCA))
            prev <- prev + qr.fitted(object$CCA$QR, object$CCA$Xbar)
    }
    if (display == "species") {
        out <- t(apply(v^2 %*% diag(eig), 1, cumsum))
        if (addprevious)
            out <- out + colSums(prev^2)
    } else {
        out <- matrix(0, nrow(u), ncol(u))
        for (i in seq_len(ncol(u))) {
            mat <- tcrossprod(u[,i], v[,i])
            out[,i] <- rowSums(mat^2) * eig[i]
        }
        out <- t(apply(out, 1, cumsum))
        if (addprevious)
            out <- out + prev
    }
    out
}

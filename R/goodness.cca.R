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
    ## Handle dbrda objects
    if (!inherits(object, "dbrda"))
        eig <- eig[eig > 0]
    if (inherits(object, "dbrda")) {
        u <- cbind(u, object[[model]][["imaginary.u"]])
        display <- "dbrda"
    }
    ## get the total variation
    All <- if (is.null(object$CCA)) object$CA$Xbar else object$CCA$Xbar
    if (!is.null(object$pCCA))
        All <- All + object$pCCA$Fit
    tot <- switch(display,
                  "species" = colSums(All^2),
                  "sites" = rowSums(All^2),
                  "dbrda" = diag(All))
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
            prev <- object$pCCA$Fit
        else
            prev <- 0
        if (model == "CA" && !is.null(object$CCA))
            prev <- prev + qr.fitted(object$CCA$QR, object$CCA$Xbar)
    }
    if (display == "species") {
        out <- t(apply(v^2 %*% diag(eig), 1, cumsum))
        if (addprevious)
            out <- out + colSums(prev^2)
        dimnames(out) <- dimnames(v)
    } else if (display == "sites") {
        out <- matrix(0, nrow(u), ncol(u))
        if (addprevious)
            mat <- prev
        else
            mat <- 0
        for (i in seq_len(ncol(u))) {
            mat <- tcrossprod(u[,i], v[,i]) * sqrt(eig[i]) + mat
            out[,i] <- rowSums(mat^2)
        }
        dimnames(out) <- dimnames(u)
    } else { # dbrda
        out <- matrix(0, nrow(u), ncol(u))
        mat <- 0
        for (i in seq_len(ncol(u))) {
            mat <- u[,i]^2 * eig[i] + mat
            out[,i] <- mat
        }
        dimnames(out) <- dimnames(u)
    }
    out/tot
}

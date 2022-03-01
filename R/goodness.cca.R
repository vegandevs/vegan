`goodness.cca` <-
    function (object, choices, display = c("species", "sites"),
              model = c("CCA", "CA"),
              summarize = FALSE, addprevious = FALSE, ...)
{
    display <- match.arg(display)
    model <- match.arg(model)
    if (!inherits(object, "cca"))
        stop("can be used only with objects inheriting from 'cca'")
    ## See stressplot.capscale for a model to implement goodness() --
    ## this can be done, but we don't care to do it for now and we
    ## just disable this.
    if (inherits(object, "capscale"))
        stop("not implemented for 'capscale'")
    if (inherits(object, c("capscale", "dbrda")) && display == "species")
        stop(gettextf("cannot analyse species with '%s'", object$method))
    v <- sqrt(weights(object, display="species")) * object[[model]]$v
    if (is.null(v))
        stop(gettextf("model = '%s' does not exist", model))
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
    All <- ordiYbar(object, "initial")
    if (is.null(All))
        stop("old style result object: update() your model")
    tot <- switch(display,
                  "species" = colSums(All^2),
                  "sites" = rowSums(All^2),
                  "dbrda" = diag(All))
    ## take only chosen axes within the component
    if (!missing(choices)) {
        choices <- choices[choices <= length(eig)]
        if (display != "dbrda")
            v <- v[, choices, drop = FALSE]
        if (display %in% c("sites", "dbrda"))
            u <- u[, choices, drop = FALSE]
        eig <- eig[choices]
    }
    if (addprevious) {
        if (!is.null(object$pCCA))
            prev <- ordiYbar(object, "pCCA")
        else
            prev <- 0
        if (model == "CA" && !is.null(object$CCA))
            prev <- prev + ordiYbar(object, "CCA")
    }
    if (display == "species") {
        if (length(eig) == 1)
            out <- v^2 * eig
        else
            out <- t(apply(v^2 %*% diag(eig, nrow=length(eig)), 1, cumsum))
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
    if (summarize)
        out <- out[, ncol(out)]
    out/tot
}

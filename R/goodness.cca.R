`goodness.cca` <-
    function (object, display = c("species", "sites"), choices,
              model = c("CCA", "CA"),
              summarize = FALSE, addprevious = FALSE, ...)
{
    display <- match.arg(display)
    model <- match.arg(model)
    if (!inherits(object, "cca"))
        stop("can be used only with objects inheriting from 'cca'")
    if (inherits(object, c("capscale", "dbrda")) && display == "species")
        stop(gettextf("cannot analyse species with '%s'", object$method))
    what <- if(display == "species") "v" else "u"
    w <- weights(object, display = display)
    pCCA <- object$pCCA$Fit
    CA <- object[[model]][[what]]
    if (is.null(CA))
        stop(gettextf("model = '%s' does not exist", model))
    eig <- object[[model]]$eig
    if (!inherits(object, "dbrda"))
        eig <- eig[eig > 0]
    ## imaginary dimensions for dbrda
    if (inherits(object, "dbrda"))
        CA <- cbind(CA, object[[model]][["imaginary.u"]])
    ## take only chosen axes within the component
    if (!missing(choices)) {
        choices <- choices[choices <= ncol(CA)]
        CA <- CA[, choices, drop = FALSE]
        eig <- eig[choices]
    }
    att <- attributes(CA)
    if (!is.null(pCCA)) {
        if (display == "sites")
            pCCA <- t(pCCA)
        if (inherits(object, "dbrda"))
            pCCA <- diag(pCCA)
        else
            pCCA <- diag(crossprod(pCCA))
    }
    CA <- t(apply(
        diag(w, length(w)) %*% CA^2 %*% diag(eig, length(eig)),
        1, cumsum))
    ## rank=1 solutions comes out transposed: back transpose
    if (length(eig) == 1)
        CA <- t(CA)
    totals <- inertcomp(object, display = display)
    comps <- colnames(totals)
    ## statistic: explained variation
    tot <- rowSums(totals)
    if (addprevious) {
        if ("pCCA" %in% comps)
            CA <- sweep(CA, 1, totals[,"pCCA"], "+")
        if (model == "CA" && "CCA" %in% comps)
            CA <- sweep(CA, 1, totals[, "CCA"], "+")
    }
    CA <- sweep(CA, 1, tot, "/")
    ## out
    attributes(CA) <- att
    if (summarize)
        CA <- CA[,ncol(CA)]
    CA
}

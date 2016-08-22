`goodness.cca` <-
    function (object, display = c("species", "sites"), choices,
              model = c("CCA", "CA"),
              statistic = c("explained", "distance"),
              summarize = FALSE, addprevious = FALSE, ...)
{
    display <- match.arg(display)
    model <- match.arg(model)
    statistic <- match.arg(statistic)
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
    if (inherits(object, "rda"))
        nr <- nobs(object) - 1
    else
        nr <- 1
    if (!is.null(pCCA)) {
        if (display == "sites")
            pCCA <- t(pCCA)
        if (inherits(object, "dbrda"))
            pCCA <- diag(pCCA)
        else
            pCCA <- diag(crossprod(pCCA))/nr
    }
    CA <- t(apply(
        diag(w, length(w)) %*% CA^2 %*% diag(eig, length(eig)),
        1, cumsum))
    ## rank=1 solutions comes out transposed: back transpose
    if (length(eig) == 1)
        CA <- t(CA)
    totals <- inertcomp(object, display = display)
    comps <- colnames(totals)
    if (statistic == "explained") {
        tot <- rowSums(totals)
        if (addprevious) {
            if ("pCCA" %in% comps)
                CA <- sweep(CA, 1, totals[,"pCCA"], "+")
            if (model == "CA" && "CCA" %in% comps)
                CA <- sweep(CA, 1, totals[, "CCA"], "+")
        }
        CA <- sweep(CA, 1, tot, "/")
    } else {
        if ("CA" %in% comps)
            tot <- totals[,"CA"]
        else
            tot <- 0
        if (model == "CCA" && "CCA" %in% comps)
            tot <- totals[,"CCA"] + tot
        CA <- sweep(-CA, 1, tot, "+")
        CA[CA < 0] <- 0
        CA <- sqrt(CA)
        CA <- sweep(CA, 1, sqrt(w), "/")
    }
    attributes(CA) <- att
    if (summarize)
        CA <- CA[,ncol(CA)]
    CA
}

`goodness.cca` <-
    function (object, display = c("species", "sites"), choices,
              model = c("CCA", "CA"), statistic = c("explained", "distance"),
              summarize = FALSE, addprevious = FALSE, ...)
{
    model <- match.arg(model)
    display <- match.arg(display)
    if (inherits(object, "capscale") && display == "species") 
        stop("display = \"species\" not available for 'capscale'")
    if (inherits(object, "rda"))
        NR <- nobs(object) - 1
    else
        NR <- 1
    if (is.null(object$CCA)) 
        model <- "CA"
    if (is.null(object[[model]]) || object[[model]]$rank == 0) 
        stop("model ", model, " is not available")
    statistic <- match.arg(statistic)
    if (inherits(object, "rda"))
        cs <- 1
    else {
        cs <-
            if (display == "species") object$colsum else object$rowsum
    }
    lambda2 <- sqrt(object[[model]]$eig)
    ## collect contributions to the variation and scores
    ptot <- ctot <- rtot <- 0
    if (display == "species") {
        if (!is.null(object$pCCA))
            ptot <- diag(crossprod(object$pCCA$Fit)) / NR
        if (!is.null(object$CCA)) {
            Xbar <- qr.fitted(object$CCA$QR, object$CCA$Xbar)
            ctot <- diag(crossprod(Xbar)) / NR
        }
        if (!is.null(object$CA))
            rtot <- diag(crossprod(object$CA$Xbar)) / NR
        v <- sweep(object[[model]]$v, 2, lambda2, "*")
    }
    else {
        if (!is.null(object$pCCA))
            ptot <- diag(tcrossprod(object$pCCA$Fit)) / NR
        if (!is.null(object$CCA)) {
            Xbar <- qr.fitted(object$CCA$QR, object$CCA$Xbar)
            ctot <- diag(tcrossprod(Xbar)) / NR
        }
        if (!is.null(object$CA))
            rtot <- diag(tcrossprod(object$CA$Xbar)) / NR
        v <- sweep(object[[model]]$u, 2, lambda2, "*")
    }
    v <- sweep(v, 1, sqrt(cs), "*")
    if (ncol(v) > 1)
        vexp <- t(apply(v^2, 1, cumsum))
    else
        vexp <- v^2
    if (!missing(choices)) 
        vexp <- vexp[, choices, drop = FALSE]
    if (statistic == "explained") {
        tot <- ptot + ctot + rtot
        if (addprevious) {
            if (!is.null(object$pCCA))
                vexp <- sweep(vexp, 1, ptot, "+")
            if (model == "CA" && !is.null(object$CCA))
                vexp <- sweep(vexp, 1, ctot, "+")
        }
        vexp <- sweep(vexp, 1, tot, "/")
    }
    else {
        tot <- rtot
        if (model == "CCA")
            tot <- tot + ctot
        vexp <- sweep(-(vexp), 1, tot, "+")
        vexp[vexp < 0] <- 0
        vexp <- sqrt(vexp)
        vexp <- sweep(vexp, 1, sqrt(cs), "/")
    }
    if (summarize) 
        vexp <- vexp[, ncol(vexp)]
    vexp
}

`goodness2.cca` <-
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
    eig <- object[[model]]$eig
    eig <- eig[eig > 0]
    ## imaginary dimensions for dbrda
    if (inherits(object, "dbrda"))
        CA <- cbind(CA, object[[model]][["imaginary.u"]])
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
    CA <- t(apply(diag(w) %*% CA^2 %*% diag(eig), 1,
                  cumsum))
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
    if (summarize)
        CA <- CA[,ncol(CA)]
    attributes(CA) <- att
    CA
}

`goodness.cca` <-
    function (object, display = c("species", "sites"), choices,
              model = c("CCA", "CA"), statistic = c("explained", "distance"),
              summarize = FALSE, addpartial = TRUE, ...) 
{
    model <- match.arg(model)
    display <- match.arg(display)
    if (is.null(object$CCA)) 
        model <- "CA"
    if (is.null(object[[model]]) || object[[model]]$rank == 0) 
        stop("model ", model, " is not available")
    statistic <- match.arg(statistic)
    cs <- if (display == "species") object$colsum else object$rowsum
    lambda2 <- sqrt(object[[model]]$eig)
    ## collect contributions to the variation and scores
    ptot <- ctot <- rtot <- 0
    if (display == "species") {
        if (!is.null(object$pCCA))
            ptot <- diag(crossprod(object$pCCA$Fit))
        if (!is.null(object$CCA)) {
            Xbar <- qr.fitted(object$CCA$QR, object$CCA$Xbar)
            ctot <- diag(crossprod(Xbar))
        }
        if (!is.null(object$CA))
            rtot <- diag(crossprod(object$CA$Xbar))
        v <- sweep(object[[model]]$v, 2, lambda2, "*")
    }
    else {
        if (!is.null(object$pCCA))
            ptot <- diag(tcrossprod(object$pCCA$Fit))
        if (!is.null(object$CCA)) {
            Xbar <- qr.fitted(object$CCA$QR, object$CCA$Xbar)
            ctot <- diag(tcrossprod(Xbar))
        }
        if (!is.null(object$CA))
            rtot <- diag(tcrossprod(object$CA$Xbar))
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
        if (addpartial && model == "CCA" && !is.null(object$pCCA))
            vexp <- sweep(vexp, 1, ptot, "+")
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

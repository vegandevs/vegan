`goodness.rda` <-
    function (object, display = c("species", "sites"), choices,
              model = c("CCA", "CA"), statistic = c("explained", "distance"),
              summarize = FALSE, ...) 
{
    model <- match.arg(model)
    display <- match.arg(display)
    if (inherits(object, "capscale") && display == "species") 
        stop("display = \"species\" not available for 'capscale'")
    if (is.null(object$CCA)) 
        model <- "CA"
    if (is.null(object[[model]]) || object[[model]]$rank == 0) 
        stop("model ", model, " is not available")
    statistic <- match.arg(statistic)
    lambda2 <- sqrt(object[[model]]$eig)
    ## collect contributions to the variation and scores
    ptot <- ctot <- rtot <- 0
    if (display == "species") {
        if (!is.null(object$pCCA))
            ptot <- diag(crossprod(object$pCCA$Fit))/(nrow(object$pCCA$Fit)-1)
        if (!is.null(object$CCA)) {
            Xbar <- qr.fitted(object$CCA$QR, object$CCA$Xbar)
            ctot <- diag(crossprod(Xbar)) / (nrow(Xbar) - 1)
        }
        if (!is.null(object$CA))
            rtot <- diag(crossprod(object$CA$Xbar)) / (nrow(object$CA$Xbar) - 1)
        v <- sweep(object[[model]]$v, 2, lambda2, "*")
    }
    else {
        if (!is.null(object$pCCA))
            ptot <- diag(tcrossprod(object$pCCA$Fit))/(nrow(object$pCCA$Fit)-1)
        if (!is.null(object$CCA)) {
            Xbar <- qr.fitted(object$CCA$QR, object$CCA$Xbar)
            ctot <- diag(tcrossprod(Xbar)) / (nrow(Xbar) - 1)
        }
        if (!is.null(object$CA))
            rtot <- diag(tcrossprod(object$CA$Xbar)) / (nrow(object$CA$Xbar) - 1)
        v <- sweep(object[[model]]$u, 2, lambda2, "*")
    }
    if (ncol(v) > 1)
        vexp <- t(apply(v^2, 1, cumsum))
    else
        vexp <- v^2
    if (!missing(choices)) 
        vexp <- vexp[, choices, drop = FALSE]
    if (statistic == "explained") {
        tot <- ptot + ctot + rtot
        vexp <- sweep(vexp, 1, tot, "/")
    }
    else {
        tot <- rtot
        if (model == "CCA")
            tot <- tot + ctot
        vexp <- sweep(-(vexp), 1, tot, "+")
        vexp[vexp < 0] <- 0
    }
    if (summarize) 
        vexp <- vexp[, ncol(vexp)]
    vexp
}


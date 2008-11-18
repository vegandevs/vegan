`predict.rda` <-
    function (object, newdata, type = c("response", "wa", "sp", "lc"), 
              rank = "full", model = c("CCA", "CA"), scaling = FALSE, ...) 
{
    type <- match.arg(type)
    model <- match.arg(model)
    if (model == "CCA" && is.null(object$CCA)) 
        model <- "CA"
    take <- object[[model]]$rank
    if (rank != "full") 
        take <- min(take, rank)
    if (is.null(object$CCA)) 
        tmp <- object$CA$Xbar
    else tmp <- object$CCA$Xbar
    cent <- attr(tmp, "scaled:center")
    scal <- attr(tmp, "scaled:scale")
    scaled.PCA <- !is.null(scal)
    nr <- nrow(tmp) - 1
    u <- object[[model]]$u[, 1:take, drop = FALSE]
    v <- object[[model]]$v[, 1:take, drop = FALSE]
    w <- object[[model]]$wa[, 1:take, drop = FALSE]
    if (is.null(w)) 
        w <- u
    slam <- diag(sqrt(object[[model]]$eig[1:take] * nr), nrow = take)
    if (type == "response") {
        if (inherits(object, "capscale")) 
            stop("Prediction of 'response' not available in capscale")
        if (!is.null(object$pCCA)) 
            warning("Conditional ('partial') component ignored")
        if (take > 0) 
            out <- u %*% slam %*% t(v)
        if (!is.null(scal)) 
            out <- sweep(out, 2, scal, "*")
        out <- sweep(out, 2, cent, "+")
    }
    else if (type == "lc") {
        if (model == "CA") 
            stop("'lc' scores not available for unconstrained ordination")
        if (!missing(newdata)) {
            if (is.null(object$terminfo)) 
                E <- as.matrix(newdata)
            else {
                d <- ordiParseFormula(formula(object), newdata, 
                                      object$terminfo$xlev)
                E <- cbind(d$Z, d$Y)
            }
            Q <- object[[model]]$QR
            p1 <- Q$pivot[1:Q$rank]
            E <- sweep(E, 2, c(object$pCCA$envcentre, object$CCA$envcentre), 
                       "-")
            u <- E[, p1, drop = FALSE] %*% coef(object)[p1, , 
                         drop = FALSE]
            u <- u[, 1:take, drop = FALSE]
        }
        out <- u
        if (scaling) {
            tot <- sqrt(object$tot.chi * nr)
            lam <- list(diag(slam)/tot, 1, sqrt(diag(slam)/tot))[[abs(scaling)]]
            out <- sqrt(tot) * sweep(out, 2, lam, "*")
        }
    }
    else if (type == "wa") {
        if (!missing(newdata)) {
            if (inherits(object, "capscale")) 
                stop("'wa' scores not available in capscale with 'newdata'")
            if (!is.null(object$pCCA)) 
                stop("No 'wa' scores available (yet) in partial RDA")
            Xbar <- as.matrix(newdata)
            Xbar <- sweep(Xbar, 2, cent, "-")
            if (!is.null(scal)) 
                Xbar <- sweep(Xbar, 2, scal, "/")
            w <- Xbar %*% v
            w <- sweep(w, 2, diag(slam), "/")
        }
        out <- w
        if (scaling) {
            tot <- sqrt(object$tot.chi * nr)
            lam <- list(diag(slam)/tot, 1, sqrt(diag(slam)/tot))[[abs(scaling)]]
            out <- sqrt(tot) * sweep(out, 2, lam, "*")
        }
    }
    else if (type == "sp") {
        if (inherits(object, "capscale")) 
            warning("'sp' scores may be meaningless in 'capscale'")
        if (!missing(newdata)) {
            Xbar <- as.matrix(newdata)
            Xbar <- scale(Xbar, center = TRUE, scale = scaled.PCA)
            if (!is.null(object$pCCA)) 
                Xbar <- qr.resid(object$pCCA$QR, Xbar)
            v <- t(Xbar) %*% u
            v <- sweep(v, 2, diag(slam), "/")
        }
        out <- v
        if (scaling) {
            tot <- sqrt(object$tot.chi * nr)
            scal <- list(1, diag(slam)/tot, sqrt(diag(slam)/tot))[[abs(scaling)]]
            out <- sqrt(tot) * sweep(out, 2, scal, "*")
        }
    }
    out
}

"predict.cca" <-
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
    rs <- object$rowsum
    cs <- object$colsum
    gtot <- object$grand.tot
    u <- object[[model]]$u[, 1:take, drop = FALSE]
    v <- object[[model]]$v[, 1:take, drop = FALSE]
    w <- object[[model]]$wa[, 1:take, drop = FALSE]
    if (is.null(w)) 
        w <- u
    slam <- diag(sqrt(object[[model]]$eig[1:take]), nrow = take)
    if (type == "response") {
        Xbar <- 0
        if (take > 0) 
            Xbar <- u %*% slam %*% t(v)
        if (!is.null(object$pCCA)) 
            warning("Conditional ('partial') component ignored")
        rc <- outer(rs, cs)
        out <- (Xbar + 1) * rc * gtot
    }
    else if (type == "lc") {
        if (model == "CA") 
            stop("'lc' scores not available for unconstrained ordination")
        if (!missing(newdata)) {
            if (is.null(object$terminfo))
                E <- as.matrix(newdata)
            else {
                d <- ordiParseFormula(formula(object), newdata, object$terminfo$xlev)
                E <- cbind(d$Z, d$Y)
            }
            E <- sweep(E, 2, c(object$pCCA$envcentre, object$CCA$envcentre), 
                       "-")
            Q <- object[[model]]$QR
            p1 <- Q$pivot[1:Q$rank]
            u <- E[, p1, drop = FALSE] %*% coef(object)[p1, , 
                         drop = FALSE]
            u <- u[, 1:take, drop = FALSE]
        }
        out <- u
        if (scaling) {
            scal <- list(diag(slam), 1, sqrt(diag(slam)))[[abs(scaling)]]
            out <- sweep(out, 2, scal, "*")
            if (scaling < 0) {
                scal <- sqrt(1/(1 - diag(slam)^2))
                out <- sweep(out, 2, scal, "*")
            }
        }
    }
    else if (type == "wa") {
        if (!missing(newdata)) {
            if (!is.null(object$pCCA)) 
                stop("No 'wa' scores available (yet) in partial CCA")
            Xbar <- as.matrix(newdata)
            exclude.spec <- attr(object[[model]]$v, "na.action")
            if (!is.null(exclude.spec)) 
                Xbar <- Xbar[, -exclude.spec]
            rs <- rowSums(Xbar)
            Xbar <- (Xbar - outer(rs, cs))/sqrt(outer(rs, cs))
            v <- sweep(v, 1, sqrt(cs), "*")
            w <- sweep(Xbar %*% v, 1, sqrt(rs), "/")
            w <- sweep(w, 2, diag(slam), "/")
        }
        out <- w
        if (scaling) {
            scal <- list(diag(slam), 1, sqrt(diag(slam)))[[abs(scaling)]]
            out <- sweep(out, 2, scal, "*")
            if (scaling < 0) {
                scal <- sqrt(1/(1 - diag(slam)^2))
                out <- sweep(out, 2, scal, "*")
            }
        }
    }
    else if (type == "sp") {
        if (!missing(newdata)) {
            Xbar <- as.matrix(newdata)
            cs <- colSums(Xbar)
            Xbar <- (Xbar - outer(rs, cs))/sqrt(outer(rs, cs))
            if (!is.null(object$pCCA)) 
                Xbar <- qr.resid(object$pCCA$QR, Xbar)
            u <- sweep(u, 1, sqrt(rs), "*")
            v <- sweep(t(Xbar) %*% u, 1, sqrt(cs), "/")
            v <- sweep(v, 2, diag(slam), "/")
        }
        out <- v
        if (scaling) {
            scal <- list(1, diag(slam), sqrt(diag(slam)))[[abs(scaling)]]
            out <- sweep(out, 2, scal, "*")
            if (scaling < 0) {
                scal <- sqrt(1/(1 - diag(slam)^2))
                out <- sweep(out, 2, scal, "*")
            }
        }
    }
    out
}

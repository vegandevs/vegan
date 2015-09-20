`predict.cca` <-
    function (object, newdata, type = c("response", "wa", "sp", "lc", "working"), 
              rank = "full", model = c("CCA", "CA"), scaling = "none",
              hill = FALSE, ...) 
{
    type <- match.arg(type)
    model <- match.arg(model)
    if (model == "CCA" && is.null(object$CCA)) 
        model <- "CA"
    take <- object[[model]]$rank
    if (take == 0)
        stop("model ", dQuote(model), " has rank 0")
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
    ## process scaling arg, scaling used later so needs to be a numeric
    scaling <- scalingType(scaling = scaling, hill = hill)
    if (type %in%  c("response", "working")) {
        Xbar <- 0
        if (!missing(newdata)) {
            if (NROW(u) == NROW(newdata))
                u <- predict(object, type = if(model == "CCA") "lc" else "wa",
                             newdata = newdata, rank = take)
            else
                warning(gettextf("'newdata' ignored: it must have the same number of rows as the original community data with type = '%s'", type))
        }
        if (take > 0) 
            Xbar <- u %*% slam %*% t(v)
        rc <- outer(rs, cs)
        if (type == "response") 
            out <- (Xbar + 1) * rc * gtot
        else                 # type == "working"
            out <- Xbar * sqrt(rc)
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
        if (scaling) {      # implicit conversion "none" == 0 == FALSE
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
            nm <- rownames(v)
            if (!is.null(nm)) { # Got rownames: keep only species with scores
                if (!all(nm %in% colnames(newdata)))
                    stop("'newdata' does not have named columns matching one or more the original columns")
                newdata <-  newdata[, nm, drop = FALSE]
            } else { #Rownames are NULL: still try to remove exclude.spec
                exclude.spec <- attr(object[[model]]$v, "na.action")
                if (!is.null(exclude.spec)) 
                    Xbar <- Xbar[, -exclude.spec]
            }
            Xbar <- as.matrix(newdata)
            rs <- rowSums(Xbar)
            Xbar <- (Xbar - outer(rs, cs))/sqrt(outer(rs, cs))
            v <- sweep(v, 1, sqrt(cs), "*")
            w <- sweep(Xbar %*% v, 1, sqrt(rs), "/")
            w <- sweep(w, 2, diag(slam), "/")
        }
        out <- w
        if (scaling) {      # implicit conversion "none" == 0 == FALSE
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
            nm <- rownames(u)
            if (!is.null(nm)) {
                if (!all(nm %in% rownames(newdata)))
                    stop("'newdata' does not have named rows matching one or more of the original rows")
                newdata <- newdata[nm, , drop = FALSE]
            }
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
        if (scaling) {      # implicit conversion "none" == 0 == FALSE
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

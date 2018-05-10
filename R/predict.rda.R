### predict.rda handles rda plus distance-based capscale and
### dbrda. Distance-based methods have some limitations:
###
### - type = "response" returns dissimilarities (and ignores imaginary dims)
### - Euclidean distances type = "working" give type = "response"
### - there are no meaningful species scores
### - WA scores with newdata cannot be calculated in capscale.
### - only type = "response", "working" and "lc" work with dbrda
### - only type = "lc" can be used with newdata with dbrda

`predict.rda` <-
    function (object, newdata, type = c("response", "wa", "sp", "lc", "working"),
              rank = "full", model = c("CCA", "CA"), scaling = "none",
              correlation = FALSE, const, ...)
{
    ## not vegan rda, but intended for klaR:::predict.rda?
    if (!("CA" %in% names(object)))
        stop(gettextf("%s is not a vegan rda object",
                      sQuote(deparse(substitute(object)))))
    type <- match.arg(type)
    model <- match.arg(model)
    if (model == "CCA" && is.null(object$CCA))
        model <- "CA"
    if (inherits(object, "dbrda"))
        take <- object[[model]]$poseig
    else
        take <- object[[model]]$rank
    if (take == 0)
        stop(gettextf("model '%s' has rank 0", model))
    if (rank != "full")
        take <- min(take, rank)
    if (!inherits(object, "dbrda")) {
        ## the ifs are only needed to cope with pre-2.5-0 vegan: now
        ## we always have Ybar, but earlier we needed to check whether
        ## we had CA or CCA Xbar
        if (!is.null(object$Ybar)) {
            cent <- attr(object$Ybar, "scaled:center")
            scal <- attr(object$Ybar, "scaled:scale")
        } else { # needed for vegan-2.4 compatibility
            if (is.null(object$CCA))
                tmp <- object$CA$Xbar
            else tmp <- object$CCA$Xbar
            cent <- attr(tmp, "scaled:center")
            scal <- attr(tmp, "scaled:scale")
        }
        scaled.PCA <- !is.null(scal)
    }
    nr <- nobs(object) - 1
    u <- object[[model]]$u[, 1:take, drop = FALSE]
    w <- object[[model]]$wa[, 1:take, drop = FALSE]
    if (is.null(w))
        w <- u
    if (!inherits(object, "dbrda")) {
        v <- object[[model]]$v[, 1:take, drop = FALSE]
    }
    ## process scaling arg, scaling used later so needs to be a numeric
    scaling <- scalingType(scaling = scaling, correlation = correlation)
    if (type %in% c("wa","sp","lc")) {
        slam <- sqrt(object[[model]]$eig[1:take]/object$tot.chi)
        if (scaling && missing(const))
            const <- sqrt(sqrt((nobs(object)-1) * object$tot.chi))
    } else {
        slam <- diag(sqrt(object[[model]]$eig[1:take]), nrow = take)
    }

    if (type %in% c("response", "working")) {
        if (!missing(newdata)) {
            u <- predict(object, type = if(model == "CCA") "lc" else "wa",
                         newdata = newdata, rank = take)
        }
        if (inherits(object, c("capscale", "dbrda"))) {
            if (take > 0) {
                out <- u %*% slam
                if (type == "response") {
                    out <- dist(out)
                    if (!is.null(object$ac)) {
                        if (object$add == "lingoes")
                            out <- sqrt(out^2 - 2 * object$ac)
                        else if (object$add == "cailliez")
                            out <- out - object$ac
                        else
                            stop("unknown euclidifying adjustment")
                    }
                    if (object$sqrt.dist)
                        out <- out^2
                }
            }
        } else {
            if (take > 0)
                out <- u %*% slam %*% t(v)
            else {
                out <- matrix(0, nrow = nrow(u), ncol = nrow(v))
                rownames(out) <- rownames(u)
                colnames(out) <- rownames(v)
            }
            if (type == "response") {
                if (!is.null(scal))
                    out <- sweep(out, 2, scal, "*")
                out <- out * sqrt(nr)
                out <- sweep(out, 2, cent, "+")
            }
        }
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
        if (scaling) {   # implicit coercion 0 == FALSE, other == TRUE
            lam <- list(slam, 1, sqrt(slam))[[abs(scaling)]]
            out <- const * sweep(out, 2, lam, "*")
        }
    }
    else if (type == "wa") {
        if (!missing(newdata)) {
            if (inherits(object, c("capscale", "dbrda")))
                stop(gettextf("'wa' scores not available in %s with 'newdata'",
                     object$method))
            if (!is.null(object$pCCA))
                stop("no 'wa' scores available (yet) in partial RDA")
            nm <- rownames(v)
            if (!is.null(nm)) {
                if (!all(nm %in% colnames(newdata)))
                    stop("'newdata' does not have named columns matching one or more the original columns")
                newdata <-  newdata[, nm, drop = FALSE]
            }
            Xbar <- as.matrix(newdata)
            Xbar <- sweep(Xbar, 2, cent, "-")
            Xbar <- Xbar / sqrt(nr)
            if (!is.null(scal)) {
                nz <- scal > 0
                Xbar[,nz] <- sweep(Xbar[,nz], 2, scal[nz], "/")
            }
            w <- Xbar %*% v
            w <- sweep(w, 2, slam, "/") / sqrt(object$tot.chi)
        }
        out <- w
        if (scaling) {   # implicit coercion 0 == FALSE, other == TRUE
            lam <- list(slam, 1, sqrt(slam))[[abs(scaling)]]
            out <- const * sweep(out, 2, lam, "*")
        }
    }
    else if (type == "sp") {
        if (inherits(object, "capscale"))
            warning("'sp' scores may be meaningless in 'capscale'")
        if (inherits(object, "dbrda"))
            stop("'sp' scores are not available in 'dbrda'")
        if (!missing(newdata)) {
            nm <- rownames(u)
            if (!is.null(nm)) {
                if (!all(nm %in% rownames(newdata)))
                    stop("'newdata' does not have named rows matching one or more of the original rows")
                newdata <- newdata[nm, , drop = FALSE]
            }
            Xbar <- as.matrix(newdata)
            Xbar <- scale(Xbar, center = TRUE, scale = scaled.PCA)
            Xbar <- Xbar/sqrt(nr)
            if (!is.null(object$pCCA))
                Xbar <- qr.resid(object$pCCA$QR, Xbar)
            v <- t(Xbar) %*% u
            v <- sweep(v, 2, slam, "/") / sqrt(object$tot.chi)
        }
        out <- v
        if (scaling) {   # implicit coercion 0 == FALSE, other == TRUE
            lam <- list(1, slam, sqrt(slam))[[abs(scaling)]]
            out <- const * sweep(out, 2, lam, "*")
            if (scaling < 0) { # correlation=TRUE
                out <- out / object$colsum
                out <- out * sqrt(object$tot.chi / (nobs(object)-1))
            }
        }
    }
    out
}

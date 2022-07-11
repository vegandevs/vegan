### prediction mainly fails with detrended (and rescaled) analysis
### beyond axis 1 or with 'newdata'. It could be possible to work out
### approximate prediction, but that would require orthogonalizing the
### scores and using additive eigenvalues (evals.ortho): see
### decorana() for orthogonalization.

`predict.decorana` <-
    function (object, newdata, type = c("response", "sites", "species"),
              rank = 4, ...)
{
    type <- match.arg(type)
    u <- object$rproj[, 1:rank, drop = FALSE]
    v <- object$cproj[, 1:rank, drop = FALSE]
    orig <- object$origin[1:rank]
    u <- sweep(u, 2, orig, "-")
    v <- sweep(v, 2, orig, "-")
    rs <- object$aidot
    cs <- object$adotj
    tot <- sum(rs)
    rs <- rs/tot
    cs <- cs/tot
    lam <- object$evals[1:rank]
    if (type == "response") {
        if (!object$ira && rank > 1)
            stop("prediction is unavailable in detrended CA beyond first axis")
        Xbar <- 0
        if (rank > 0) {
            if (!object$ira) {
                tmp <- sweep(v, 1, sqrt(cs), "*")
                tmp <- rbind(tmp, sweep(u, 1, sqrt(rs), "*"))
                rot <- svd(tmp)$v
                v <- v %*% rot
                u <- u %*% rot
                fac <- colSums(sweep(v^2, 1, cs, "*"))
                lam <- (fac - 1)/fac
            }
            Xbar <- u %*% diag(1 - lam, nrow = rank) %*% t(v)
        }
        out <- (Xbar + 1) * outer(rs, cs) * tot
    }
    else if (type == "sites") {
        if (!missing(newdata)) {
            Xbar <- as.matrix(newdata)
            if (!is.null(object$v))
                Xbar <- sweep(Xbar, 2, object$v, "*")
            rs <- rowSums(Xbar)
            Xbar <- (Xbar - outer(rs, cs))/sqrt(outer(rs, cs))
            v <- sweep(v, 1, sqrt(cs), "*")
            u <- sweep(Xbar %*% v, 1, sqrt(rs), "/")
        }
        out <- u
    }
    else if (type == "species") {
        if (!missing(newdata)) {
            if (!object$ira && rank > 1)
                stop("type = 'species' not available in detrended CA with 'newdata'")
            if (object$iresc)
                stop("type = 'species' not available in rescaled DCA with 'newdata'")
            Xbar <- as.matrix(newdata)
            cs <- colSums(Xbar)
            Xbar <- (Xbar - outer(rs, cs))/sqrt(outer(rs, cs))
            u <- sweep(u, 1, sqrt(rs), "*")
            v <- sweep(t(Xbar) %*% u, 1, sqrt(cs), "/")
            v <- sweep(v, 2, lam, "/")
        }
        out <- v
    }
    out
}

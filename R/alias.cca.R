`alias.cca` <-
    function (object, names.only = FALSE, ...) 
{
    if (names.only)
        return(object$CCA$alias)
    CompPatt <- function(x, ...) {
        x[abs(x) < 1e-06] <- 0
        class(x) <- "mtable"
        x[abs(x) < 1e-06] <- NA
        x
    }
    Model <- object$terms
    attributes(Model) <- NULL
    value <- list(Model = Model)
    R <- object$CCA$QR$qr
    R <- R[1:min(dim(R)), , drop = FALSE]
    R[lower.tri(R)] <- 0
    d <- dim(R)
    rank <- object$CCA$QR$rank
    p <- d[2]
    value$Complete <- if (is.null(p) || rank == p) 
        NULL
    else {
        p1 <- 1:rank
        X <- R[p1, p1]
        Y <- R[p1, -p1, drop = FALSE]
        beta12 <- as.matrix(qr.coef(qr(X), Y))
        CompPatt(t(beta12))
    }
    class(value) <- "listof"
    value
}


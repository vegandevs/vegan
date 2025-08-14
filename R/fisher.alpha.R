`fisher.alpha` <-
    function (x, MARGIN = 1, ...)
{
    x <- as.matrix(x, rownames.force = TRUE)
    if(ncol(x) == 1)
        x <- t(x)
    sol <- apply(x, MARGIN, fisherfit)
    out <-  unlist(lapply(sol, function(x) x$estimate))
    if (FALSE) {
        out <- list(alpha = out)
        out$se <- unlist(lapply(sol, function(x) sqrt(diag(solve(x$hessian)))[1]))
        out$df.residual <- unlist(lapply(sol, df.residual))
        out$code <- unlist(lapply(sol, function(x) x$code))
        out <- as.data.frame(out)
    }
    out
}

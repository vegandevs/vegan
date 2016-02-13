`raupcrick` <-
    function(comm, null = "r1", nsimul = 999, chase = FALSE, ...)
{
    comm <- as.matrix(comm)
    comm <- ifelse(comm > 0, 1, 0)
    ## 'tri' is a faster alternative to as.dist(): it takes the lower
    ## diagonal, but does not set attributes of a "dist" object
    N <- nrow(comm)
    tri <- matrix(FALSE, N, N)
    tri <- row(tri) > col(tri)
    ## function(x) designdist(x, "J", terms="binary") does the same,
    ## but is much slower
    sol <- oecosimu(comm, function(x) tcrossprod(x)[tri], method = null,
                    nsimul = nsimul,
                    alternative = if (chase) "less" else "greater",
                    ...)
    ## Chase et al. way, or the standard way
    if (chase)
        out <- 1 - sol$oecosimu$pval
    else
        out <- sol$oecosimu$pval
    ## set attributes of a "dist" object
    attributes(out) <- list("class"=c("raupcrick", "dist"), "Size"=N,
                            "Labels" = rownames(comm), "call" =
                            match.call(), "Diag" = FALSE, "Upper" = FALSE,
                            "method" = "raupcrick")
    out
}

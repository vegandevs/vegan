`spantree` <-
    function (d, toolong = 0)
{
    d <- as.dist(d)
    n <- attr(d, "Size")
    labels <- labels(d)
    dis <- .C(primtree, dist = as.double(d), toolong = as.double(toolong),
              n = as.integer(n), val = double(n + 1),
              dad = integer(n + 1), NAOK = TRUE)
    out <- list(kid = dis$dad[2:n] + 1, dist = dis$val[2:n],
                labels = labels, n = n, call = match.call())
    class(out) <- "spantree"
    out
}

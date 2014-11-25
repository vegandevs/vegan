`spantree` <-
    function (d, toolong = 0) 
{
    dis <- as.dist(d)
    n <- attr(dis, "Size")
    labels <- labels(dis)
    dis <- .C("primtree", dist = as.double(dis), toolong = as.double(toolong), 
              n = as.integer(n), val = double(n + 1),
              dad = integer(n + 1), NAOK = TRUE, PACKAGE = "vegan")
    out <- list(kid = dis$dad[2:n] + 1, dist = dis$val[2:n],
                labels = labels, n = n, call = match.call())
    class(out) <- "spantree"
    out
}

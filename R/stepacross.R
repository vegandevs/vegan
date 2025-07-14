`stepacross` <-
    function (dis, path = "shortest", toolong = 1, trace = TRUE, ...)
{
    path <- match.arg(path, c("shortest", "extended"))
    if (!inherits(dis, "dist"))
        dis <- as.dist(dis)
    oldatt <- attributes(dis)
    n <- attr(dis, "Size")
    if (path == "shortest")
        dis <- .C(dykstrapath, dist = as.double(dis), n = as.integer(n),
                  as.double(toolong), as.integer(trace),
                  out = double(length(dis)), NAOK = TRUE, PACKAGE = "vegan")$out
    else dis <- .C(C_stepacross, dis = as.double(dis), as.integer(n),
                   as.double(toolong), as.integer(trace), NAOK = TRUE,
                   PACKAGE = "vegan")$dis
    if("maxdist" %in% oldatt)
        oldatt$maxdist <- NA
    attributes(dis) <- oldatt
    attr(dis, "method") <- paste(attr(dis, "method"), path)
    dis
}

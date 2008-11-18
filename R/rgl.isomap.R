`rgl.isomap` <-
    function(x, web = "white", ...)
{
    require(rgl) || stop("requires package 'rgl'")
    ordirgl(x, ...)
    z <- scores(x, ...)
    net <- x$net
    for (i in 1:nrow(net))
        rgl.lines(z[net[i,],1], z[net[i,],2], z[net[i,],3], color=web)
}

## Extract the points in the hull as a one matrix
`scores.ordihull` <-
    function(x, ...)
{
    out <- NULL
    for(i in 1:length(x))
        out <- rbind(out, x[[i]])
    hulls <- rep(names(x), sapply(x, function(z) NROW(z)))
    attr(out, "hulls") <- hulls
    out
}

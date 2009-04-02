`plot.meandist` <-
    function(x, cluster = "average",  ...) 
{
    n <- attr(x, "n")
    cl <- hclust(as.dist(x), method = cluster, members = n)
    cl <- as.dendrogram(cl, hang = 0)
    w <- diag(x)[labels(cl)]
    tr <- unlist(dendrapply(cl, function(n) attr(n, "height")))
    root <- attr(cl, "height")
    plot(cl, ylim = range(c(w, tr, root), na.rm = TRUE), leaflab = "none", ...)
    for (i in 1:length(w)) segments(i, tr[i], i, w[i])
    pos <- ifelse(w < tr, 1, 3)
    pos[is.na(pos)] <- 1
    w[is.na(w)] <- tr[is.na(w)]
    text(1:length(w), w, labels = labels(cl), pos = pos, srt = 0)
}


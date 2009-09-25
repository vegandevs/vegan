`plot.meandist` <-
    function(x, kind = c("dendrogram", "histogram"),  cluster = "average", ylim,
             axes = TRUE, ...) 
{
    kind <- match.arg(kind)
    n <- attr(x, "n")
    if (kind == "dendrogram") {
        cl <- hclust(as.dist(x), method = cluster, members = n)
        cl <- as.dendrogram(cl, hang = 0)
        w <- diag(x)[labels(cl)]
        tr <- unlist(dendrapply(cl, function(n) attr(n, "height")))
        root <- attr(cl, "height")
        if (missing(ylim))
            ylim <- range(c(w, tr, root), na.rm = TRUE)
        plot(cl, ylim = ylim, leaflab = "none", axes = axes, ...)
        for (i in 1:length(w)) segments(i, tr[i], i, w[i])
        pos <- ifelse(w < tr, 1, 3)
        pos[is.na(pos)] <- 1
        w[is.na(w)] <- tr[is.na(w)]
        text(1:length(w), w, labels = labels(cl), pos = pos, srt = 0)
    } else {
        w <- diag(x)
        tr <- rep(summary(x)$B, length(w))
        if (missing(ylim))
            ylim <- range(c(w, tr), na.rm = TRUE)
        plot(1:length(w), tr, ylim = ylim, axes = FALSE, xlab = "", ylab = "",
             type = "l", ...)
        if (axes)
            axis(2, ...)
        for (i in 1:length(w)) segments(i, tr, i, w[i])
        pos <- ifelse(w < tr, 1, 3)
        pos[is.na(pos)] <- 1
        text(1:length(w), w, labels = names(n), pos = pos, srt = 0,
             xpd = TRUE, ...)
    }
}


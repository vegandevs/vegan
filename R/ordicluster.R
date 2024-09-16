`ordicluster` <-
    function (ord, cluster, prune = 0, display = "sites",
              w, col = 1, draw = c("segments", "none"), ...)
{
    if(missing(w))
        w <- if(is.atomic(ord)) attr(ord, "weights")
             else weights(ord, display = display)
    mrg <- cluster$merge
    ord <- scores(ord, display = display, ...)
    if (ncol(ord) > 2)
        ord <- ord[, 1:2, drop = FALSE]
    if (nrow(mrg) != nrow(ord) - 1)
        stop("dimensions do not match in 'ord' and 'cluster'")
    if (is.null(w) || length(w) == 1) w <- rep(w, nrow(ord))
    n <- if (is.null(w)) rep(1, nrow(ord)) else w
    noden <- numeric(nrow(mrg) - prune)
    go <- matrix(0, nrow(mrg) - prune, 2)
    ## recycle colours for points and prepare to get node colours
    col <- rep(col, length = nrow(ord))
    col <- col2rgb(col)/255
    nodecol <- matrix(NA, nrow(mrg) - prune, 3)
    seg.coords <- matrix(NA, nrow = nrow(mrg) - prune, ncol = 4)
    for (i in seq_len(nrow(mrg) - prune)) {
        a <- mrg[i,1]
        b <- mrg[i,2]
        one <- if (a < 0) ord[-a,] else go[a,]
        two <- if (b < 0) ord[-b,] else go[b,]
        n1 <- if (a < 0) n[-a] else noden[a]
        n2 <- if (b < 0) n[-b] else noden[b]
        xm <- weighted.mean(c(one[1],two[1]), w = c(n1,n2))
        ym <- weighted.mean(c(one[2],two[2]), w = c(n1,n2))
        go[i,] <- c(xm,ym)
        noden[i] <- n1 + n2
        colone <- if (a < 0) col[,-a] else nodecol[a,]
        coltwo <- if (b < 0) col[,-b] else nodecol[b,]
        nodecol[i,] <- (n1 * colone + n2 * coltwo)/noden[i]
        seg.coords[i, ] <- c(one[1], one[2], two[1], two[2])
    }
    colnames(seg.coords) <- c("x1","y1","x2","y2")
    ## are we plotting?
    draw <- match.arg(draw)
    if (isTRUE(all.equal(draw, "segments"))) {
        ordiArgAbsorber(seg.coords[,1L], seg.coords[,2L],
                        seg.coords[,3L], seg.coords[,4L],
                        col = rgb(nodecol),
                        FUN = segments, ...)
    }
    colnames(go) <- c("x","y")
    seg.coords <- cbind(as.data.frame(seg.coords), col = rgb(nodecol))
    out <- structure(list(scores = cbind(go, "w" = noden),
                          segments = seg.coords), class = "ordicluster")
    invisible(out)
}

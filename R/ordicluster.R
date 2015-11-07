`ordicluster` <-
    function (ord, cluster, prune=0, display="sites",
              w = weights(ord, display), col = 1, ...)
{
    weights.default <- function(object, ...) NULL
    w <- eval(w)
    mrg <- cluster$merge
    ord <- scores(ord, display = display, ...)
    if (nrow(mrg) != nrow(ord) - 1)
        stop("Dimensions do not match in 'ord' and 'cluster'")
    if (length(w) == 1) w <- rep(w, nrow(ord))
    n <- if (is.null(w)) rep(1, nrow(ord)) else w
    noden <- numeric(nrow(ord))
    go <- ord
    ## recycle colours for points and prepare to get node colours
    col <- rep(col, length = nrow(ord))
    col <- col2rgb(col)/255
    nodecol <- matrix(NA, nrow(mrg) - prune, 3)
    for (i in 1: (nrow(mrg) - prune)) {
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
        ordiArgAbsorber(one[1], one[2], two[1], two[2],
                        col = rgb(t(nodecol[i,])),
                        FUN = segments, ...)
    }
    invisible(cbind(go, "w"=noden))
}

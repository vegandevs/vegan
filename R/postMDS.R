"postMDS" <-
    function (X, dist, pc = TRUE, center = TRUE, halfchange = TRUE, 
              threshold = 0.8, nthreshold = 10, plot = FALSE, ...) 
{
    Size <- attributes(dist)$Size
    if (any(attributes(X)$names == "points")) 
        x <- X$points
    else x <- as.matrix(X)
    if (center) 
        x <- scale(x, scale = FALSE)
    if (pc) {
        dn <- dimnames(x)
        x <- prcomp(x, center = center)$x
        dimnames(x) <- dn
    }
    if (halfchange) {
        dist <- as.vector(dist)
        ordi <- as.vector(vegdist(x, "euclidean"))
        take <- dist < threshold
        if (sum(take) < nthreshold) {
            warning("skipping half-change scaling: too few points below threshold")
            halfchange <- FALSE
        }
        else {
            k <- coef(lm(dist[take] ~ ordi[take]))
            names(k) <- NULL
            hc <- (1 - k[1])/2/k[2]
            x <- x/hc
        }
    }
    if (plot && halfchange) {
        cross.lim <- 45
        if (Size > cross.lim) 
            pch <- "."
        else pch <- "+"
        orange <- range(c(ordi, 0, 1))
        drange <- range(c(dist, 0, 1))
        plot(orange, drange, type = "n", xlab = "Ordination distance", 
             ylab = "Community dissimilarity")
        points(ordi[take], dist[take], pch = pch, col = "blue")
        points(ordi[!take], dist[!take], pch = pch, col = "gray")
        abline(h = threshold)
        abline(h = k[1])
        hclevel <- (1 - k[1])/2 + k[1]
        segments(0, hclevel, hc, hclevel, col = "red", lwd = 2)
        arrows(hc, hclevel, hc, 0, col = "red", lwd = 2)
        arrows(0, k[1], 0, hclevel, col = "red", code = 3)
        arrows(0, hclevel, 0, 1, col = "red", code = 3)
        j <- 0.02
        text(0 + j, threshold + j, "Threshold", adj = c(0, 0))
        text(0 + j, k[1] + j, "Replicate dissimilarity", adj = c(0, 
                                                         0))
        text(0 + j, hclevel + j, "Half-change", adj = c(0, 0))
        abline(k, col = "blue", lwd = 2)
    }
    attr(x, "centre") <- center
    attr(x, "pc") <- pc
    attr(x, "halfchange") <- halfchange
    if (any(attributes(X)$names == "points")) 
        X$points <- x
    else X <- x
    X
}

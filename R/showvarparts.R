`showvarparts` <-
    function(parts, labels, bg = NULL, alpha=63, Xnames, id.size=1.2, ...)
{
    rad <- 0.725
    ## Default names
    if (missing(Xnames))
        Xnames <- paste("X", seq_len(parts), sep="")
    ## transparent fill colours
    if (!is.null(bg)) {
        bg <- rgb(t(col2rgb(bg)), alpha = alpha, maxColorValue = 255)
        if (length(bg) < parts)
            bg <- rep(bg, length.out = parts)
    }
    ## centroids of circles (parts < 4) or individual fractions (parts
    ## == 4)
    cp <- switch(parts,
                 matrix(c(0,0), ncol=2, byrow=TRUE),
                 matrix(c(0,0, 1,0), ncol=2, byrow=TRUE),
                 matrix(c(0,0, 1,0, 0.5, -sqrt(3/4)), ncol=2, byrow=TRUE),
                 structure(
                     c(-1.2, -0.6, 0.6, 1.2, -0.7, 0, -0.7, 0, 0.7, 0.7,
                       0.3, -0.4, 0.4, -0.3, 0, 0, 0.7, 0.7, 0, 0.3, 0.4,
                       -0.6,-1.2, -0.6, 0.3, -0.7, 0, 0, -0.7, -0.4),
                     .Dim = c(15L, 2L))
                 )
    ## positions of labels for cp (circles, ellipses)
    pos.names <-
        switch(parts,
               matrix(c(0,0), ncol=2, byrow=TRUE),
               matrix(c(0,1,0.8,0.8),2,2),
               matrix(c(-0.4,1.4,0.05,0.65,0.65,-1.5), 3, 2),
               matrix(c(-1.25,-0.85,0.85,1.25,0.6,1.05,1.05,0.6),4,2)
               )
    pos <- switch(parts,
                  0,
                  c(2, 4),
                  c(2, 4, 2),
                  c(2, 2, 4, 4))
    ## plot limits
    if (parts < 4) {
        xlim <- range(cp[,1]) + c(-rad, rad)
        ylim <- range(cp[,2]) + c(-rad, rad)
    } else {
        xlim <- c(-1.7, 1.7)
        ylim <- c(-1.7, 1.1)
    }
    ## plot
    plot(cp, axes=FALSE, xlab="", ylab="", asp=1, type="n",
         xlim = xlim, ylim = ylim)
    ## See if text fits the space
    xspace <- par("usr")[1:2]
    sw <- strwidth(Xnames, cex = id.size) + strwidth("m", cex = id.size)
    pt <- pos.names[,1]
    pt <- ifelse(pos == 2, pt - sw, pt + sw)
    if (any(pt < xspace[1]) || any(pt > xspace[2])) {
        message("labels do not fit the space, switching to X1...")
        Xnames <- paste0("X", seq_len(parts))
    }
    text(pos.names, labels=Xnames[seq_len(parts)], cex=id.size, pos = pos)
    box()
    if (parts < 4) {
        symbols(cp, circles = rep(rad, min(parts,3)), inches = FALSE,
                add=TRUE, bg = bg, ...)
    } else {
        ## Draw ellipses with veganCovEllipse. Supply 2x2
        ## matrix(c(d,a,a,d), 2, 2) which defines an ellipse of
        ## semi-major axis length sqrt(d+a) semi-minor axis sqrt(d-a).
        d <- 1
        a <- 1/sqrt(2)
        ## Small ellipses X2, X3 at the centroid
        e2 <- veganCovEllipse(matrix(c(d,-a,-a,d), 2, 2))
        e3 <- veganCovEllipse(matrix(c(d, a, a,d), 2, 2))
        ## wider ellipses X1, X4 at sides going through the centroid
        L <- d+a
        W <- (sqrt(L) - sqrt(d-a))^2
        d <- (L+W)/2
        a <- (L-W)/2
        cnt <- sqrt(W/2)
        e1 <- veganCovEllipse(matrix(c(d,-a,-a,d), 2, 2), c(-cnt, -cnt))
        e4 <- veganCovEllipse(matrix(c(d, a, a,d), 2, 2), c( cnt, -cnt))
        polygon(rbind(e1,NA,e2,NA,e3,NA,e4), col = bg, ...)
    }

    ## label fractions
    nlabs <- switch(parts, 2, 4, 8, 16)
    if (missing(labels))
        labels <- paste("[", letters[1:nlabs], "]", sep="")
    if (length(labels) != nlabs)
        stop(gettextf("needs %d labels, but input has %d",
                      nlabs, length(labels)))
    switch(parts,
           text(0,0, labels[-nlabs], ...),
           text(rbind(cp, colMeans(cp)), labels[-nlabs], ...),
           text(rbind(cp, colMeans(cp[1:2,]), colMeans(cp[2:3,]),
                      colMeans(cp[c(1,3),]), colMeans(cp)), labels[-nlabs], ...),
           text(cp, labels[-nlabs], ...)
           )
    xy <- par("usr")
    text(xy[2] - 0.05*diff(xy[1:2]), xy[3] + 0.05*diff(xy[3:4]),
         paste("Residuals =", labels[nlabs]), pos = 2, ...)
    invisible()
}

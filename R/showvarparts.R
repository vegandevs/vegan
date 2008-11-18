"showvarparts" <-
function(parts = 2, labels, ...)
{
    rad <- 0.725
    cp <- switch(parts,
                 c(0,0),
                 c(0,0, 1,0),
                 c(0,0, 1,0, 0.5, -sqrt(3/4)),
                 c(-0.5,0.3, 0.5, 0.3, 0, -sqrt(3/4)+0.3)
                 )
    cp <- matrix(cp, ncol=2, byrow=TRUE)
    plot(cp, axes=FALSE, xlab="", ylab="", asp=1, type="n", 
         xlim = (range(cp[,1]) + c(-rad, rad)),
         ylim = (range(cp[,2]) + c(-rad, rad)))
    box()
    symbols(cp, circles = rep(rad, min(parts,3)), inches = FALSE, add=TRUE, ...)
    if (parts == 4) {
        symbols(0, 0.2, rectangles=cbind(1, 0.5), inch=FALSE, add=TRUE, ...)
        symbols(sqrt(1/2), -sqrt(3/4)+0.2, rectangles=cbind(0.5,0.3),
                inch=FALSE, add=TRUE, ...)
    }
    nlabs <- switch(parts, 2, 4, 8, 16)
    if (missing(labels))
        labels <- paste("[", letters[1:nlabs], "]", sep="")
    if (length(labels) != nlabs)
        stop("needs ", nlabs, " labels, but input has" , length(labels))
    switch(parts,
           text(0,0, labels[-nlabs], ...),
           text(rbind(cp[1,], colMeans(cp), cp[2,]), labels[-nlabs], ...),
           text(rbind(cp, colMeans(cp[1:2,]), colMeans(cp[2:3,]),
                      colMeans(cp[c(1,3),]), colMeans(cp)), labels[-nlabs], ...),
           text(rbind(1.4*cp, c(0.8, -sqrt(3/4)+0.2),
                      colMeans(cp[1:2,]) + c(0,0.25),
                      colMeans(cp[2:3,]), colMeans(cp[c(1,3),]),
                      cp[1,] + c(0.1,0), cp[2,] -c(0.1,0),
                      c(0.6, -sqrt(3/4)+0.2), colMeans(cp[1:2,]),
                      colMeans(cp)-c(0,0.12), colMeans(cp[2:3,]) + c(0,0.14),
                      colMeans(cp[c(1,3),]) + c(0, 0.14),
                      colMeans(cp) + c(0,0.08)),
                labels[-nlabs], ...)
           )
    xy <- par("usr")
    text(xy[2] - 0.05*diff(xy[1:2]), xy[3] + 0.05*diff(xy[3:4]),
         paste("Residuals =", labels[nlabs]), pos = 2, ...)
    invisible()
}


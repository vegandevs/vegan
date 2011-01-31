`rgl.renyiaccum` <-
    function(x, rgl.height = 0.2,  ...)
{
    require(rgl) || stop("requires packages 'rgl'")
    y <- x[,,1] * rgl.height
    rgl.min = 0
    rgl.max = max(y)
    xp <- seq(0, 1, len = nrow(y))
    z <- seq(0, 1, len = ncol(y))
    ylim <- 1000 * range(y)
    ylen <- ylim[2] - ylim[1] + 1
    colorlut <- rainbow(ylen)
    col <- colorlut[1000*y-ylim[1]+1]
    rgl.bg(color = "white")
    rgl.surface(xp, z, y, color=col)
    y <- x[,,5] * rgl.height
    ##rgl.surface(xp,z,y,color="grey", alpha=0.3)
    rgl.surface(xp, z, y,  color="black", front="lines", back="lines")
    y <- x[,,6] * rgl.height
    ##rgl.surface(xp,z,y,color="grey",alpha=0.3)
    rgl.surface(xp, z, y, color="black", front="lines", back="lines")
    y <- x[,,6]*0 + rgl.min
    rgl.surface(xp, z, y, alpha=0)
    y <- x[,,6] * 0 + rgl.max
    rgl.surface(xp, z, y, alpha=0)
    labs <- pretty(c(rgl.min, range(x)))
    rgl.bbox(color="#333377", emission="#333377", specular="#3333FF", shininess=5, alpha=0.8,
             zlen=0, xlen=0, yat = rgl.height*labs, ylab=labs) 
    rgl.texts(0, rgl.min, 0.5, "Scale", col = "darkblue")
    rgl.texts(0.5, rgl.min, 0, "Sites", col="darkblue")
}

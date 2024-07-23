### Modelled after maptools:::pointLabel.
`ordipointlabel` <-
    function(x, display = c("sites", "species"), choices = c(1,2), col=c(1,2),
             pch=c("o","+"), font = c(1,1), cex=c(0.8, 0.8), add = FALSE,
             select, ...)
{
    xy <- list()
    ## Some 'scores' accept only one 'display': a workaround
    for (nm in display)
        xy[[nm]] <- scores(x, display = nm, choices = choices, ...)
    ##xy <- scores(x, display = display, choices = choices, ...)
    ## remove `select`ed observations from scores as per text.cca
    ## only useful if we are displaying only one set of scores
    if(!missing(select)) {
        if(isTRUE(all.equal(length(display), 1L))) {
            xy[[1]] <- .checkSelect(select, xy[[1]])
        } else {
            warning("'select' does not apply when plotting more than one set of scores--\n'select' was ignored")
        }
    }
    if (length(display) > 1) {
        col <- rep(col, sapply(xy, nrow))
        pch <- rep(pch, sapply(xy, nrow))
        font <- rep(font, sapply(xy, nrow))
        cex <- rep(cex, sapply(xy, nrow))
        tmp <- xy[[1]]
        for (i in 2:length(display))
            tmp <- rbind(tmp, xy[[i]])
        xy <- tmp
    }
    else {
        xy <- xy[[1]]
        if (length(col) < nrow(xy))
            col <- rep(col[1], nrow(xy))
        if (length(pch) < nrow(xy))
            pch <- rep(pch[1], nrow(xy))
        if (length(cex) < nrow(xy))
            cex <- rep(cex[1], length = nrow(xy))
        if (length(font) < nrow(xy))
            font <- rep(font[1], length = nrow(xy))
    }
    if (!add)
        pl <- ordiplot(xy, display = "sites", type="n")
    labels <- rownames(xy)
    em <- strwidth("m", cex = min(cex), font = min(font))
    ex <- strheight("x", cex = min(cex), font = min(font))
    ltr <- em*ex
    ## bounding box: strwidth/height do not accept vector cex and font
    ## and we loop
    box <- matrix(0, nrow(xy), 2)
    for (i in seq_len(nrow(xy))) {
        box[i,1] <- strwidth(labels[i], cex = cex[i], font = font[i]) +
            strwidth("m", cex = cex[i], font = font[i])
        box[i,2] <- strheight(labels[i], cex = cex[i], font = font[i]) +
            strheight("x", cex = cex[i], font = font[i])
    }
    ## offset: 1 up, 2..4 sides, 5..8 corners
    makeoff <- function(pos, lab) {
        cbind(c(0,1,0,-1,0.9,0.9,-0.9,-0.9)[pos] * lab[,1]/2,
              c(1,0,-1,0,0.8,-0.8,-0.8,0.8)[pos] * lab[,2]/2)
    }
    ## amount of overlap
    overlap <- function(xy1, off1, xy2, off2) {
        pmax(0, pmin(xy1[,1] + off1[,1]/2, xy2[,1] + off2[,1]/2)
             -pmax(xy1[,1] - off1[,1]/2, xy2[,1] - off2[,1]/2)) *
              pmax(0, pmin(xy1[,2] + off1[,2]/2, xy2[,2] + off2[,2]/2)
             -pmax(xy1[,2] - off1[,2]/2, xy2[,2] - off2[,2]/2))
    }
    ## indices of overlaps in lower triangular matrix
    n <- nrow(xy)
    j <- as.vector(as.dist(row(matrix(0, n, n))))
    k <- as.vector(as.dist(col(matrix(0, n, n))))
    ## Find labels that may overlap...
    maylap <- overlap(xy[j,], 2*box[j,], xy[k,], 2*box[k,]) > 0
    ## ... and work only with those
    j <- j[maylap]
    k <- k[maylap]
    jk <- sort(unique(c(j,k)))
    ## SANN: no. of iterations & starting positions
    nit <- min(48 * length(jk), 10000)
    pos <- rep(1, n)
    ## Criterion: overlap + penalty for positions other than directly
    ## above and especially for corners
    fn <- function(pos) {
        off <- makeoff(pos, box)
        val <- sum(overlap(xy[j,]+off[j,], box[j,], xy[k,]+off[k,], box[k,]))
        val <- val/ltr + sum(pos>1)*0.1 + sum(pos>4)*0.1
    }
    ## Move a label of one point
    gr <- function(pos) {
        take <- sample(jk, 1)
        pos[take] <- sample((1:8)[-pos[take]], 1)
        pos
    }
    ## Simulated annealing
    sol <- optim(par = pos, fn = fn, gr = gr, method="SANN",
                 control=list(maxit=nit))
    if (!add)
        ##points(xy, pch = pch, col = col, cex=cex, ...)
        ordiArgAbsorber(xy, pch = pch, col = col, cex = cex, FUN = points,
                        ...)
    lab <- xy + makeoff(sol$par, box)
    ##text(lab, labels=labels, col = col, cex = cex, font = font,  ...)
    ordiArgAbsorber(lab, labels=labels, col = col, cex = cex, font = font,
                    FUN = text, ...)
    pl <- list(points = xy)
    pl$labels <- lab
    attr(pl$labels, "font") <- font
    args <- list(tcex = cex, tcol = col, pch = pch, pcol = col,
                 pbg = NA, pcex = cex)
    pl$args <- args
    pl$par <- par(no.readonly = TRUE)
    pl$dim <- par("din")
    attr(pl, "optim") <- sol
    class(pl) <- c("ordipointlabel", "orditkplot", class(pl))
    invisible(pl)
}

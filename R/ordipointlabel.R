### Modelled after maptools:::pointLabel.
`ordipointlabel` <-
    function(x, display = c("sites", "species"), choices = c(1,2), col=c(1,2),
             pch=c("o","+"), font = c(1,1), cex=c(0.7, 0.7),
             add = inherits(x, "ordiplot"), labels, bg, select, ...)
{
    xy <- list()
    ## Some 'scores' accept only one 'display': a workaround
    for (nm in display)
        xy[[nm]] <- scores(x, display = nm, choices = choices, ...)

    if (length(display) > 1) {
        ld <- length(display)
        col <- rep(rep(col, length=ld), sapply(xy, nrow))
        pch <- rep(rep(pch, length=ld), sapply(xy, nrow))
        font <- rep(rep(font, length=ld), sapply(xy, nrow))
        cex <- rep(rep(cex, length=ld), sapply(xy, nrow))
        if (!missing(bg))
            fill <- rep(rep(bg, length=ld), sapply(xy, nrow))
        xy <- do.call(rbind, xy)
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
        if (!missing(bg) && length(bg)  < nrow(xy))
            fill <- rep(bg[1], length = nrow(xy))
    }
    if (!add)
        pl <- ordiplot(xy, display = "sites", type="n")

    ## keep only 'select': applicable only if there was only one
    ## set of scores
    if(!missing(select)) {
        if(identical(length(display), 1L)) {
            xy <- .checkSelect(select, xy)
        } else
            warning("'select' can be used only with one set of scores: ignoring 'select'")
    }
    if (!missing(labels)) {
            if (length(labels) != nrow(xy)) {
            stop(gettextf(
                "you need  %d labels but arg 'labels' only had %d: arg ignored",
                nrow(xy), length(labels)))
        }
        rownames(xy) <- labels
    } else {
        labels <- rownames(xy)
    }
    ## our algorithm needs at least three items
    if (nrow(xy) < 3) {
        warning(gettextf("optimization needs at least three items, you got %d",
                         nrow(xy)))
        return(invisible(x))
    }
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
    nit <- min(64 * length(jk), 10000)
    pos <- ifelse(xy[,2] > 0, 1, 3)
    ## Criterion: overlap + penalty for moving towards origin and also
    ## for corners. Penalty is mild: max 1 ltr and one-character
    ## string > 3*ltr due to padding (em, ex) of the bounding box.
    fn <- function(pos) {
        move <- makeoff(pos, matrix(1, 1, 2))
        off <- makeoff(pos, box)
        val <- sum(overlap(xy[j,,drop=FALSE]+off[j,,drop=FALSE],
                           box[j,,drop=FALSE],
                           xy[k,,drop=FALSE]+off[k,,drop=FALSE],
                           box[k,,drop=FALSE]))
        val <- val/ltr + sum(move[,1] * xy[,1] < 0) * 0.4 +
            sum(move[,2] * xy[,2] < 0) * 0.4 +
            sum(pos > 4) * 0.2
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
    lab <- xy + makeoff(sol$par, box)
    dev.hold()
    ## draw optional lab background first so it does not cover points
    if (!missing(bg)) {
        for(i in seq_len(nrow(lab))) {
            polygon(lab[i,1] + c(-1,1,1,-1)*box[i,1]/2.2,
                    lab[i,2] + c(-1,-1,1,1)*box[i,2]/2.2,
                    col = fill[i], border = col[i], xpd = TRUE)
            ordiArgAbsorber(lab[i,1], lab[i,2], labels = labels[i],
                            col = col[i], cex = cex[i], font = font[i],
                            FUN = text, ...)
        }
    } else {
        ordiArgAbsorber(lab, labels=labels, col = col, cex = cex,
                        font = font, FUN = text, ...)
    }

    ## always plot points (heck, the function is ordi*point*label)
    ordiArgAbsorber(xy, pch = pch, col = col, cex = cex, FUN = points,
                        ...)
    ##text(lab, labels=labels, col = col, cex = cex, font = font,  ...)
    dev.flush()
    if (!inherits(x, "ordiplot"))
        pl <- list(points = xy)
    else
        pl <- x
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

### Extract labels: useful if arg labels= is given in ordipointlabel call

`labels.ordipointlabel` <-
    function(object, ...)
{
    rownames(object$labels)
}

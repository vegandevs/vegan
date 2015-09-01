`ordiellipse` <-
    function (ord, groups, display = "sites", kind = c("sd", "se"),
              conf, draw = c("lines", "polygon", "none"),
              w = weights(ord, display), col = NULL, alpha = 127,
              show.groups, label = FALSE, border=NULL,lty=NULL, lwd=NULL, ...)
{
    weights.default <- function(object, ...) NULL
    kind <- match.arg(kind)
    draw <- match.arg(draw)
    pts <- scores(ord, display = display, ...)
    ## ordiellipse only works with 2D data (2 columns)
    pts <- as.matrix(pts)
    if (ncol(pts) > 2)
        pts <- pts[ , 1:2, drop = FALSE]
    if (ncol(pts) < 2)
        stop("ordiellipse needs two dimensions")
    w <- eval(w)
    if (length(w) == 1)
        w <- rep(1, nrow(pts))
    if (is.null(w))
        w <- rep(1, nrow(pts))
    ## make semitransparent fill
    if (draw == "polygon" && !is.null(col))
        col <- rgb(t(col2rgb(col)), alpha = alpha, maxColorValue = 255)
    if (!missing(show.groups)) {
        take <- groups %in% show.groups
        pts <- pts[take, , drop = FALSE]
        groups <- groups[take]
        w <- w[take]
    }
    out <- seq(along = groups)
    inds <- names(table(groups))
    
    # fill in graphical vectors with default values if unspecified and recycles shorter vectors 
    for(arg in c("col","border","lty","lwd")){
      tmp <- mget(arg,ifnotfound=list(NULL))[[1]]
      if(is.null(tmp)) tmp <- 1
      if(length(inds) != length(tmp)) {tmp <- rep_len(tmp, length(inds))}
      assign(arg, tmp)
      
    }
    
    res <- list()
    if (label) {
        cntrs <- matrix(NA, nrow=length(inds), ncol=2)
        rownames(cntrs) <- inds
    }
    ## Remove NA scores
    kk <- complete.cases(pts) & !is.na(groups)
    for (is in inds) {
        gr <- out[groups == is & kk]
        if (length(gr)) {
            X <- pts[gr, , drop = FALSE]
            W <- w[gr]
            mat <- cov.wt(X, W)
            if (mat$n.obs == 1)
                mat$cov[] <- 0
            if (kind == "se")
                mat$cov <- mat$cov/mat$n.obs
            if (missing(conf))
                t <- 1
            else t <- sqrt(qchisq(conf, 2))
            if (mat$n.obs > 1)
                xy <- veganCovEllipse(mat$cov, mat$center, t)
            else
                xy <- X
            if (draw == "lines")
                ordiArgAbsorber(xy, FUN = lines,
                                col = if (is.null(col)) 
                                  par("fg")
                                else col[match(is, inds)],
                                lty=lty[match(is,inds)],lwd=lwd[match(is,inds)], ...)
                      
            else if (draw == "polygon") 
                ordiArgAbsorber(xy[, 1], xy[, 2], col = col[match(is, inds)], border=border[match(is,inds)],
                                lty=lty[match(is,inds)],lwd=lwd[match(is,inds)],
                                FUN = polygon,
                                ...)
            if (label && draw != "none") {
                cntrs[is,] <- mat$center
            }
            mat$scale <- t
            res[[is]] <- mat
        }
    }
    if (label && draw != "none") {
        if (draw == "lines")
            ordiArgAbsorber(cntrs[,1], cntrs[,2], labels = rownames(cntrs),
                            col = col,  FUN = text, ...)
        else 
            ordiArgAbsorber(cntrs, col = NULL,
                            FUN = ordilabel, ...)
    }
    class(res) <- "ordiellipse"
    invisible(res)
}

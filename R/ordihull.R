`ordihull` <-
    function (ord, groups, display = "sites",
              draw = c("lines", "polygon", "none"),
              col = NULL, alpha = 127, show.groups, label = FALSE,
              border=NULL, lty=NULL, lwd=NULL, ...)
      
{
    draw <- match.arg(draw)
    ## Internal function to find the polygon centre
    polycentre <- function(x) {
        n <- nrow(x)
        if (n < 4) 
            return(colMeans(x[-n, , drop = FALSE]))
        xy <- x[-n, 1] * x[-1, 2] - x[-1, 1] * x[-n, 2]
        A <- sum(xy)/2
        xc <- sum((x[-n, 1] + x[-1, 1]) * xy)/A/6
        yc <- sum((x[-n, 2] + x[-1, 2]) * xy)/A/6
        c(xc, yc)
    }
    ## Make semitransparent fill colour
    if (draw == "polygon" && !is.null(col))
        col <- rgb(t(col2rgb(col)), alpha = alpha, maxColorValue = 255)
    pts <- scores(ord, display = display, ...)
    if (!missing(show.groups)) {
        take <- groups %in% show.groups
        pts <- pts[take, , drop = FALSE]
        groups <- groups[take]
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
            X <- pts[gr,, drop = FALSE]
            hpts <- chull(X)
            hpts <- c(hpts, hpts[1])
            if (draw == "lines") 
              ordiArgAbsorber(X[hpts, ], FUN = lines, col = if (is.null(col)) 
                par("fg")
                else col[match(is, inds)], lty=lty[match(is,inds)],lwd=lwd[match(is,inds)], ...)
            else if (draw == "polygon") 
              ordiArgAbsorber(X[hpts, ], border = border[match(is,inds)],FUN = polygon, col = col[match(is, inds)], ...)
            if (label && draw != "none") {
                cntrs[is,] <- polycentre(X[hpts,])
            }
            res[[is]] <- X[hpts,]
        }
    }
    if (label && draw != "none") {
      if (draw == "lines") 
        ordiArgAbsorber(cntrs[, 1], cntrs[, 2], labels = names, 
                        col = col[match(is, inds)], FUN = text, ...)
      else ordiArgAbsorber(cntrs, labels = names, col = NULL, 
                           FUN = ordilabel, ...)
    }
    class(res) <- "ordihull"
    invisible(res)
}

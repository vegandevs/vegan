`ordihull` <-
    function (ord, groups, display = "sites",
              draw = c("lines", "polygon", "none"),
              col = NULL, alpha = 127, show.groups, label = FALSE,
              border = NULL, lty = NULL, lwd = NULL, ...)
      
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
    ## Make semitransparent fill colour; alpha should be integer
    ## 0..255, but we also handle real values < 1
    if (alpha < 1)
        alpha <- round(alpha * 255)
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
    
    ## fill in graphical vectors with default values if unspecified
    ## and recycles shorter vectors
    col.new <- border.new <- lty.new <- lwd.new <- NULL
    for(arg in c("col","border","lty","lwd")){
      tmp <- mget(arg,ifnotfound=list(NULL))[[1]]
      if(is.null(tmp))
          tmp <- ifelse(suppressWarnings(is.null(par(arg))),
                        par("fg"), par(arg))
      if(length(inds) != length(tmp))
          tmp <- rep_len(tmp, length(inds))
      assign(paste(arg,".new", sep=""), tmp)
    }
    ## default colour for "polygon" fill is "transparent", for lines
    ## is par("fg")
    if(is.null(col) && draw=="polygon")
        col.new <- rep_len("transparent", length(inds))
    else if(is.null(col) && draw=="lines")
        col.new <- rep_len(par("fg"), length(inds))
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
                ordiArgAbsorber(X[hpts, ], FUN = lines,
                                col = if (is.null(col)) 
                                          par("fg")
                                      else
                                          col.new[match(is, inds)],
                                lty = lty.new[match(is,inds)],
                                lwd = lwd.new[match(is,inds)], ...)
            else if (draw == "polygon") 

                ordiArgAbsorber(X[hpts, ],
                                border= border.new[match(is,inds)],
                                FUN = polygon,
                                col = col.new[match(is, inds)],
                                lty = lty.new[match(is,inds)],
                                lwd=lwd.new[match(is,inds)], ...)

            if (label && draw != "none") {
                cntrs[is,] <- polycentre(X[hpts,])
            }
            res[[is]] <- X[hpts,]
        }
    }
    if (label && draw != "none") {
      if (draw == "lines") 
          ordiArgAbsorber(cntrs[, 1], cntrs[, 2],
                          labels = rownames(cntrs), 
                          col = col.new[match(is, inds)],
                          FUN = text, ...)
      else ordiArgAbsorber(cntrs, labels = rownames(cntrs),
                           col = NULL, 
                           FUN = ordilabel, ...)
    }
    class(res) <- "ordihull"
    invisible(res)
}

### The function displays (ordered) heatmaps of community data. It
### copies vegemite() for handling 'use', 'sp.ind', 'site.ind' and
### 'select', but then switches to heatmap() to display the
### data. Unlike heatmap(), it does not insist on showing dendrograms,
### but only uses these for sites, and only if given as 'use'.

`tabasco` <-
    function (x, use, sp.ind = NULL, site.ind = NULL,  
              select, Rowv = TRUE, Colv = TRUE, ...) 
{
    if (any(x < 0))
        stop("function cannot be used with negative data values")
    pltree <- sptree <- NA
    if (!missing(use)) {
        if (!is.list(use) && is.vector(use)) {
            if (is.null(site.ind)) 
                site.ind <- order(use)
            if (is.null(sp.ind)) 
                sp.ind <- order(wascores(use, x))
        }
        else if (inherits(use, c("dendrogram", "hclust", "twins"))) {
            ## "twins" and "dendrogram" are treated as "dendrogram",
            ## but "hclust" is kept as "hclust": they differ in
            ## reorder()
            if (inherits(use, "twins")) {
                use <- as.dendrogram(use)
            }
            if (!is.null(site.ind))
                stop("'site.ind' cannot be used with dendrogram")
            ## The tree/dendrogam and input data must be ordered
            ## identically. It could be regarded as a "user error" if
            ## they are not, but this could be really frustrating and
            ## give obscure errors, and therefore we take care of
            ## identical ordering here
            if (inherits(use, "hclust") && !is.null(use$labels))
                x <- x[use$labels,]
            else # dendrogram
                x <- x[labels(use),]
            ## Reorder tree if Rowv specified
            if (isTRUE(Rowv)) {
                ## order by first CA axis -- decorana() is fastest
                tmp <- decorana(x, ira = 1)
                ## reorder() command is equal to all, but "dendrogram"
                ## will use unweighted mean and "hclust" weighted
                ## mean.
                use <- reorder(use, scores(tmp, dis="sites", choices = 1),
                               agglo.FUN = "mean")
            } else if (length(Rowv) > 1) {
                ## Rowv is a vector
                if (length(Rowv) != nrow(x))
                    stop(gettextf("Rowv has length %d, but 'x' has %d rows",
                                  length(Rowv), nrow(x)))
                use <- reorder(use, Rowv, agglo.FUN = "mean")
            }
            if (inherits(use, "dendrogram")) { 
                site.ind <- seq_len(nrow(x))
                names(site.ind) <- rownames(x)
                site.ind <- site.ind[labels(use)]
            } else {
                site.ind <- use$order
            }
            if (is.null(sp.ind)) 
                sp.ind <- order(wascores(order(site.ind), x))
            pltree <- use
            ## heatmap needs a "dendrogram"
            if(!inherits(pltree, "dendrogram"))
                pltree <- as.dendrogram(pltree)
        }
        else if (is.list(use)) {
            tmp <- scores(use, choices = 1, display = "sites")
            if (is.null(site.ind)) 
                site.ind <- order(tmp)
            if (is.null(sp.ind)) 
                sp.ind <- try(order(scores(use, choices = 1, 
                                           display = "species")))
            if (inherits(sp.ind, "try-error")) 
                sp.ind <- order(wascores(tmp, x))
        }
        else if (is.matrix(use)) {
            tmp <- scores(use, choices = 1, display = "sites")
            if (is.null(site.ind)) 
                site.ind <- order(tmp)
            if (is.null(sp.ind)) 
                sp.ind <- order(wascores(tmp, x))
        }
    }
    ## see if sp.ind is a dendrogram or hclust tree
    if (inherits(sp.ind, c("hclust", "dendrogram", "twins"))) {
        if (inherits(sp.ind, "twins")) {
            sp.ind <- as.dendrogram(sp.ind)
        }
        sptree <- sp.ind
        ## Reorder data to match order in the dendrogam (see 'use' above)
        if (inherits(sptree, "hclust"))
            x <- x[, sptree$labels]
        else # dendrogram
            x <- x[, labels(sptree)]
        ## Consider reordering species tree
        if (isTRUE(Colv) && !is.null(site.ind)) {
            sptree <- reorder(sptree, wascores(order(site.ind), x),
                                  agglo.FUN = "mean")
        } else if (length(Colv) > 1) {
            if (length(Colv) != ncol(x))
                stop(gettextf("Colv has length %d, but 'x' has %d columns",
                              length(Colv), ncol(x)))
            sptree <- reorder(sptree, Colv, agglo.FUN = "mean")
        }
        if (inherits(sptree, "dendrogram")) {
            sp.ind <- seq_len(ncol(x))
            names(sp.ind) <- colnames(x)
            sp.ind <- sp.ind[labels(sptree)]
        } else {
            sp.ind <- sptree$order
        }
        if (!inherits(sptree, "dendrogram"))
            sptree <- as.dendrogram(sptree)
        ## reverse: origin in the upper left corner
        sptree <- rev(sptree)
    }
    if (!is.null(sp.ind) && is.logical(sp.ind))
        sp.ind <- (1:ncol(x))[sp.ind]
    if (!is.null(site.ind) && is.logical(site.ind))
        site.ind <- (1:nrow(x))[site.ind]
    if (is.null(sp.ind)) 
        sp.ind <- 1:ncol(x)
    if (is.null(site.ind)) 
        site.ind <- 1:nrow(x)
    if (!missing(select)) {
        if (!is.na(pltree))
            stop("sites cannot be 'select'ed with dendrograms or hclust trees")
        if (!is.logical(select))
            select <- sort(site.ind) %in% select
        stake <- colSums(x[select, , drop = FALSE]) > 0
        site.ind <- site.ind[select[site.ind]]
        site.ind <- site.ind[!is.na(site.ind)]
    }
    else {
        stake <- colSums(x[site.ind, ]) > 0
    }
    sp.ind <- sp.ind[stake[sp.ind]]
    ## heatmap will reorder items by dendrogram so that we need to
    ## give indices in the unsorted order if rows or columns have a
    ## dendrogram
    if (is.na(pltree[1]))
        rind <- site.ind
    else
        rind <- sort(site.ind)
    if (is.na(sptree[1]))
        ## reverse: origin in the upper left corner
        cind <- rev(sp.ind)
    else
        cind <- sort(sp.ind)
    ## we assume t() changes data.frame to a matrix
    x <- t(x[rind, cind])
    sp.nam <- rownames(x)
    sp.len <- max(nchar(sp.nam))
    heatmap((max(x) - x), Rowv = sptree, Colv = pltree,
             scale = "none", ...)
    out <- list(sites = site.ind, species = sp.ind)
    invisible(out)
}

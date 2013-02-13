### The function displays (ordered) heatmaps of community data. It
### copies vegemite() for handling 'use', 'sp.ind', 'site.ind' and
### 'select', but then switches to heatmap() to display the
### data. Unlike heatmap(), it does not insist on showing dendrograms,
### but only uses these for sites, and only if given as 'use'.

`tabasco` <-
    function (x, use, sp.ind = NULL, site.ind = NULL,  
              select, ...) 
{
    Rowv <- Colv <- NA
    if (!missing(use)) {
        if (!is.list(use) && is.vector(use)) {
            if (is.null(site.ind)) 
                site.ind <- order(use)
            if (is.null(sp.ind)) 
                sp.ind <- order(wascores(use, x))
        }
        else if (inherits(use, "hclust")) {
            if (!is.null(site.ind))
                stop("'hclust' tree cannot be 'use'd with 'site.ind'")
            site.ind <- use$order
            if (is.null(sp.ind)) 
                sp.ind <- order(wascores(order(site.ind), x))
            Colv <- as.dendrogram(use)
        }
        else if (inherits(use, "dendrogram")) {
            if (!is.null(site.ind))
                stop("'dendrogram' cannot be 'use'd with 'site.ind'")
            site.ind <- seq_len(nrow(x))
            names(site.ind) <- rownames(x)
            site.ind <- site.ind[labels(use)]
            if (is.null(sp.ind)) 
                sp.ind <- order(wascores(order(site.ind), x))
            Colv <- use
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
    if (inherits(sp.ind, c("hclust", "dendrogram"))) {
        if (!inherits(sp.ind, "dendrogram"))
            sp.ind <- as.dendrogram(sp.ind)
        Rowv <- sp.ind
        sp.ind <- seq_len(ncol(x))
        names(sp.ind) <- colnames(x)
        sp.ind <- sp.ind[labels(Rowv)]
        ## reverse: origin in the upper left corner
        Rowv <- rev(Rowv)
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
        if (!is.na(Colv))
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
    if (is.na(Colv[1]))
        rind <- site.ind
    else
        rind <- sort(site.ind)
    if (is.na(Rowv[1]))
        ## reverse: origin in the upper left corner
        cind <- rev(sp.ind)
    else
        cind <- sort(sp.ind)
    ## we assume t() changes data.frame to a matrix
    x <- t(x[rind, cind])
    sp.nam <- rownames(x)
    sp.len <- max(nchar(sp.nam))
    heatmap((max(x) - x), Rowv, Colv,  scale = "none", ...)
    out <- list(sites = site.ind, species = sp.ind)
    invisible(out)
}

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
            site.ind <- seq_len(nrow(x))
            if (is.null(sp.ind)) 
                sp.ind <- order(wascores(order(use$order), x))
            Colv <- as.dendrogram(use)
        }
        else if (inherits(use, "dendrogram")) {
            if (!is.null(site.ind))
                stop("'dendrogram' cannot be 'use'd with 'site.ind'")
            site.ind <- seq_len(nrow(x))
            o <- seq_len(nrow(x))
            names(o) <- rownames(x)
            o <- o[labels(use)]
            if (is.null(sp.ind)) 
                sp.ind <- order(wascores(order(o), x))
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
        if (inherits(use, c("hclust", "dendrogram")))
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
    x <- x[site.ind, sp.ind]
    x <- as.matrix(x)
    x <- t(x)
    sp.nam <- rownames(x)
    sp.len <- max(nchar(sp.nam))
    heatmap((max(x) - x), Rowv, Colv,  scale = "none", ...)
    out <- list(sites = site.ind, species = sp.ind)
    invisible(out)
}

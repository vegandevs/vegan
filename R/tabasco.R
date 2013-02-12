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
            if (is.null(site.ind)) 
                site.ind <- use$order
            if (is.null(sp.ind)) 
                sp.ind <- order(wascores(order(site.ind), x))
            Colv <- as.dendrogram(use)
        }
        else if (inherits(use, "dendrogram")) {
            if (is.null(site.ind)) {
                site.ind <- 1:nrow(x)
                names(site.ind) <- rownames(x)
                site.ind <- site.ind[labels(use)]
            }
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
    if (!is.null(sp.ind) && is.logical(sp.ind))
        sp.ind <- (1:ncol(x))[sp.ind]
    if (!is.null(site.ind) && is.logical(site.ind))
        site.ind <- (1:nrow(x))[site.ind]
    if (is.null(sp.ind)) 
        sp.ind <- 1:ncol(x)
    if (is.null(site.ind)) 
        site.ind <- 1:nrow(x)
    if (!missing(select)) {
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

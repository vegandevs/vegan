`vegemite` <-
    function (x, use, scale, sp.ind = NULL, site.ind = NULL, zero = ".",
              select, diagonalize = FALSE, ...)
{
    if (!missing(use)) {
        ## derived index should be based on transformed & tabulated data
        xprime <- if (missing(scale)) x
                  else coverscale(x, scale = scale, character = FALSE)
        if (!is.list(use) && is.vector(use)) {
            if (is.null(site.ind))
                site.ind <- order(use)
            if (is.null(sp.ind))
                sp.ind <- order(wascores(use, xprime))
        }
        else if (inherits(use, c("hclust", "twins"))) {
            if (inherits(use, "twins")) {
                use <- as.hclust(use)
            }
            if (diagonalize) {
                wts <- scores(cca(xprime), choices=1, display = "wa")
                use <- reorder(use, wts)
            }
            if (is.null(site.ind))
                site.ind <- use$order
            if (is.null(sp.ind))
                sp.ind <- order(wascores(order(site.ind), xprime))
        }
        else if (inherits(use, "dendrogram")) {
            if (diagonalize) {
                wts <- scores(cca(xprime), choices=1, display = "wa")
                use <- reorder(use, wts)
            }
            if (is.null(site.ind)) {
                site.ind <- seq_len(nrow(x))
                names(site.ind) <- rownames(x)
                site.ind <- site.ind[labels(use)]
            }
            if (is.null(sp.ind))
                sp.ind <- order(wascores(order(site.ind), xprime))
        }
        else if (is.list(use)) {
            tmp <- scores(use, choices = 1, display = "sites")
            if (is.null(site.ind))
                site.ind <- order(tmp)
            if (is.null(sp.ind))
                sp.ind <- try(order(scores(use, choices = 1,
                                           display = "species")))
            if (inherits(sp.ind, "try-error"))
                sp.ind <- order(wascores(tmp, xprime))
        }
        else if (is.matrix(use)) {
            tmp <- scores(use, choices = 1, display = "sites")
            if (is.null(site.ind))
                site.ind <- order(tmp)
            if (is.null(sp.ind))
                sp.ind <- order(wascores(tmp, xprime))
        }
        else if (is.factor(use)) {
            tmp <- as.numeric(use)
            if (diagonalize) {
                ord <- scores(cca(xprime, use), choices = 1,
                              display = c("lc","wa","sp"))
                if (cor(tmp, ord$constraints, method = "spearman") < 0) {
                    ord$constraints <- -ord$constraints
                    ord$sites <- -ord$sites
                    ord$species <- -ord$species
                }
                ## order factors and sites within factor levels
                site.ind <- order(round(ord$constraints, 6), ord$sites)
                if (is.null(sp.ind))
                    sp.ind <- order(ord$species)
            }
            if (is.null(site.ind))
                site.ind <- order(tmp)
            if (is.null(sp.ind))
                sp.ind <- order(wascores(tmp, xprime))
        }
    } # end of handling 'use'
    if (!is.null(sp.ind) && is.logical(sp.ind))
        sp.ind <- seq_len(ncol(x))[sp.ind]
    if (!is.null(site.ind) && is.logical(site.ind))
        site.ind <- seq_len(nrow(x))[site.ind]
    if (is.null(sp.ind))
        sp.ind <- seq_len(ncol(x))
    if (is.null(site.ind))
        site.ind <- seq_len(nrow(x))
    if (!missing(select)) {
        if (!is.logical(select))
            select <- sort(site.ind) %in% select
        stake <- colSums(x[select, , drop = FALSE]) > 0
        site.ind <- site.ind[select[site.ind]]
        site.ind <- site.ind[!is.na(site.ind)]
    }
    else {
        stake <- colSums(x[site.ind,, drop = FALSE ]) > 0
    }
    sp.ind <- sp.ind[stake[sp.ind]]
    x <- x[site.ind, sp.ind, drop = FALSE]
    if (!missing(scale))
        x <- coverscale(x, scale, ...)
    usedscale <- attr(x, "scale")
    if (any(apply(x, 1, nchar) > 1))
        stop("cowardly refusing to use longer than one-character symbols:\nUse scale")
    x <- as.matrix(x)
    x <- t(x)
    sp.nam <- rownames(x)
    sp.len <- max(nchar(sp.nam))
    nst <- ncol(x)
    nlen <- max(nchar(colnames(x)))
    page.width <- getOption("width")
    per.page <- page.width - sp.len - 3
    istart <- seq(1, nst, by = per.page)
    iend <- pmin(istart + per.page - 1, nst)
    for (st in seq_along(istart)) {
        tbl <- apply(x[, istart[st]:iend[st], drop = FALSE],
                     1, paste, sep = "", collapse = "")
        names(tbl) <- NULL
        tbl <- gsub("0", zero, tbl)
        tbl <- cbind(sp.nam, tbl)
        st.nam <- colnames(x)[istart[st]:iend[st]]
        mathead <- matrix(" ", nrow = length(st.nam), ncol = nlen)
        for (i in seq_along(st.nam)) {
            tmp <- unlist(strsplit(st.nam[i], NULL))
            start <- nlen - length(tmp) + 1
            mathead[i, start:nlen] <- tmp
        }
        head <- cbind(apply(mathead, 2, paste, sep = "", collapse = ""))
        tbl <- rbind(cbind(matrix(" ", nrow = nrow(head), 1),
                           head), tbl)
        d <- list()
        l <- 0
        for (i in dim(tbl)) {
            d[[l <- l + 1]] <- rep("", i)
        }
        dimnames(tbl) <- d
        print(noquote(tbl))
        ## collect all pages for output table
        if (exists(".tabout", inherits = FALSE))
            .tabout[,2] <- paste0(.tabout[,2], tbl[,2])
        else
            .tabout <- tbl
    }
    out <- list(sites = site.ind, species = sp.ind, table = .tabout)
    cat(length(out$sites), "sites,", length(out$species), "species\n")
    if (!is.null(usedscale))
        cat("scale: ",  usedscale, "\n")
    invisible(out)
}

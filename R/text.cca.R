`text.cca` <-
    function (x, display = "sites", labels, choices = c(1, 2),
              scaling = "species", arrow.mul, head.arrow = 0.05, select,
              const, correlation = FALSE, hill = FALSE, ...)
{
    if (length(display) > 1)
        stop("only one 'display' item can be added in one command")
    pts <- scores(x, choices = choices, display = display, scaling = scaling,
                  const, correlation = correlation, hill = hill, tidy=FALSE,
                  droplist = FALSE)
    class(pts) <- "ordiplot"
    if (!missing(labels))
        rownames(pts[[1]]) <- labels
    if (!missing(select))
        pts[[1]] <- .checkSelect(select, pts[[1]])
    text.ordiplot(pts, what = names(pts), ...)
    invisible()
}

### utility function to extract labels used in CCA/RDA/dbRDA plots:
### you may need this if you want to set your own labels=.

`labels.cca` <-
    function(object, display, ...)
{
    if (is.null(object$CCA))
        CCA <- "CA"
    else
        CCA <- "CCA"
    switch(display,
           "sp" =,
           "species" = rownames(object[[CCA]]$v),
           "wa" =,
           "sites" =,
           "lc" = rownames(object[[CCA]]$u),
           "reg" = colnames(object[[CCA]]$QR$qr),
           "bp" = rownames(object[[CCA]]$biplot),
           "cn" = {cn <- rownames(object[[CCA]]$centroids)
                   bp <- rownames(object[[CCA]]$biplot)
                   c(cn, bp[!(bp %in% cn)]) }
           )
}

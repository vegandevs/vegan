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
    if (!missing(select))
        pts[[1]] <- .checkSelect(select, pts[[1]])
    if (!missing(labels))
        rownames(pts[[1]]) <- labels
    text.ordiplot(pts, what = names(pts), ...)
    invisible()
}

### utility function to extract labels used in CCA/RDA/dbRDA plots:
### you may need this if you want to set your own labels=.

`labels.cca` <-
    function(object, display, ...)
{
    display <- match.arg(display,
               c("sp", "species", "wa", "sites", "lc", "reg", "bp", "cn"))
    if (is.null(object$CCA) || object$CCA$rank == 0)
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

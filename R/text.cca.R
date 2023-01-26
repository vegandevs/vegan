`text.cca` <-
    function (x, display = "sites", labels, choices = c(1, 2), scaling = "species",
              arrow.mul, head.arrow = 0.05, select, const, axis.bp = FALSE,
              correlation = FALSE, hill = FALSE, ...)
{
    if (length(display) > 1)
        stop("only one 'display' item can be added in one command")
    pts <- scores(x, choices = choices, display = display, scaling = scaling,
                  const, correlation = correlation, hill = hill, tidy=FALSE)
    ## store rownames of pts for use later, otherwise if user supplies
    ## labels, the checks in "cn" branch fail and "bp" branch will
    ## be entered even if there should be no "bp" plotting
    cnam <- rownames(pts)
    if (missing(labels))
        labels <- labels.cca(x, display)
    if (!missing(select)) {
        pts <- .checkSelect(select, pts)
        labels <- labels[select]
    }
    ## centroids ("cn") have special treatment: also plot biplot
    ## arrows ("bp") for continuous variables and ordered factors.
    if (display == "cn") {
        if (!is.null(nrow(pts))) { # has "cn"
            cnlabs <- seq_len(nrow(pts))
            text(pts, labels = labels[cnlabs], ...)
        } else {
            cnlabs <- NULL
        }
        pts <- scores(x, choices = choices, display = "bp", scaling = scaling,
                      const, correlation = correlation, hill = hill,
                      tidy=FALSE)
        bnam <- rownames(pts)
        pts <- pts[!(bnam %in% cnam), , drop = FALSE]
        if (nrow(pts) == 0)
            return(invisible())
        else {
            display <- "bp"
            if (!is.null(cnlabs))
                labels <- labels[-cnlabs]
        }
    }
    ## draw arrows before adding labels
    if (display %in% c("bp", "reg", "re", "r")) {
        if (missing(arrow.mul)) {
            arrow.mul <- ordiArrowMul(pts)
        }
        pts <- pts * arrow.mul
        arrows(0, 0, pts[, 1], pts[, 2], length = head.arrow,
               ...)
        pts <- ordiArrowTextXY(pts, labels, rescale = FALSE, ...)
        if (axis.bp) {
            axis(side = 3, at = c(-arrow.mul, 0, arrow.mul),
                 labels = rep("", 3))
            axis(side = 4, at = c(-arrow.mul, 0, arrow.mul),
                 labels = c(-1, 0, 1))
        }
    }
    text(pts, labels = labels, ...)
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

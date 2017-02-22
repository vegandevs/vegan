#' Extract internal working matrices from constrained ordination object
#'
#' Function extract internal working matrices from the constrained
#' ordination object (class \code{"cca"}) in \pkg{vegan}. The function
#' returns only one model matrix type, but there is no overlap between
#' elements of typical request.
#'
#' @param x ordination object
#' @param model model to be extracted
#'
`ordiYbar` <-
    function(x, model = c("CCA", "CA", "pCCA", "partial", "initial"))
{
    model <- match.arg(model)
    Ybar <- x$Ybar
    if (model == "initial")
        return(Ybar)
    ## return NULL for missing elements
    if (model != "partial")
        if(is.null(x[[model]]))
            NULL

    ## edit Ybar -- not yet dbrda
    switch(model,
    "pCCA" = {
        Ybar <- qr.fitted(x$pCCA$QR, Ybar)
    },
    "partial" = {
        if (!is.null(x$pCCA))
            Ybar <- qr.resid(x$pCCA$QR, Ybar)
    },
    "CCA" = {
        if (!is.null(x$pCCA))
            Ybar <- qr.resid(x$pCCA$QR, Ybar)
        Ybar <- qr.fitted(x$CCA$QR, Ybar)
    },
    "CA" = {
        if (!is.null(x$CCA))
            Ybar <- qr.resid(x$CCA$QR, Ybar)
        else if (!is.null(x$pCCA))
            Ybar <- qr.resid(x$pCCA$QR, Ybar)
    })
    Ybar
}

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
    isDB <- inherits(x, "dbrda")
    if (!is.null(x$Ybar))
        Ybar <- x$Ybar
    else ## vegan < 2.5 from 2016 did not have Ybar
        stop(gettextf("outdated object from 2016 or before: use update(%s)",
                      noquote(deparse(substitute(x)))))
    if (model == "initial")
        return(Ybar)
    ## return NULL for missing elements
    if (model != "partial")
        if(is.null(x[[model]]))
            return(NULL)

    ## edit Ybar -- not yet dbrda
    switch(model,
    "pCCA" = {
        Ybar <- qr.fitted(x$pCCA$QR, Ybar)
        if (isDB)
            Ybar <- qr.fitted(x$pCCA$QR, t(Ybar))
    },
    "partial" = {
        if (!is.null(x$pCCA)) {
            Ybar <- qr.resid(x$pCCA$QR, Ybar)
            if (isDB)
                Ybar <- qr.resid(x$pCCA$QR, t(Ybar))
        }
    },
    "CCA" = {
        if (!is.null(x$pCCA)) {
            Ybar <- qr.resid(x$pCCA$QR, Ybar)
            if (isDB)
                Ybar <- qr.resid(x$pCCA$QR, t(Ybar))
        }
        Ybar <- qr.fitted(x$CCA$QR, Ybar)
        if (isDB)
            Ybar <- qr.fitted(x$CCA$QR, t(Ybar))
    },
    "CA" = {
        if (!is.null(x$CCA)) {
            Ybar <- qr.resid(x$CCA$QR, Ybar)
            if (isDB)
                Ybar <- qr.resid(x$CCA$QR, t(Ybar))
        }
        else if (!is.null(x$pCCA)) {
            Ybar <- qr.resid(x$pCCA$QR, Ybar)
            if (isDB)
                Ybar <- qr.resid(x$pCCA$QR, t(Ybar))
        }
    })
    Ybar
}


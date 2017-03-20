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
    else
        return(vegan24Xbar(x, model))
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

#' Extract internal working matrices of cca objects for vegan 2.4

`vegan24Xbar` <-
    function(x, model = c("CCA", "CA", "pCCA", "partial", "initial"))
{
    model <- match.arg(model)
    ## return NULL for missing elements
    if (model %in% c("CCA", "CA", "pCCA"))
        if(is.null(x[[model]]))
            return(NULL)
    Ybar <- NULL
    ## initial working data is not saved and must be reconstructed in
    ## partial models, but this cannot be done in all cases.
    if (model == "initial") {
        if (inherits(x, "dbrda")) {
            ## NULL in partial models
            if (is.null(x$pCCA)) {
                Ybar <- x$CCA$Xbar
                if (is.null(Ybar))
                    Ybar <- x$CA$Xbar
            }
        } else {
            ## this will ignore imaginary components in capscale
            if (is.null(x$CCA))
                Ybar <- x$CA$Xbar
            else
                Ybar <- x$CCA$Xbar
            if (!is.null(x$pCCA))
                Ybar <- Ybar + x$pCCA$Fit
        }
        return(Ybar)
    }
    ## several components are already stored in the result object and
    ## we just fetch those (only CCA needs work)
    switch(model,
    "pCCA" =
        Ybar <- x$pCCA$Fit,
    "partial" = {
        Ybar <- x$CCA$Xbar
        if (is.null(Ybar))
            Ybar <- x$CA$Xbar
    },
    "CCA" = {
        Ybar <- qr.fitted(x$CCA$QR, x$CCA$Xbar)
        if (inherits(x, "dbrda"))
            Ybar <- qr.fitted(x$CCA$QR, t(Ybar))
    },
    "CA" =
        Ybar <- x$CA$Xbar
    )
    Ybar
}

#' Add Species Scores to Ordination Results
#'
#' @param object Ordination object
#' @param value Community data
#'

`sppscores<-` <-
    function(object, value)
{
    UseMethod("sppscores<-")
}

## dbrda

`sppscores<-.dbrda` <-
    function(object, value)
{
    object$vdata <- deparse(substitute(value))
    value <- scale(value, center = TRUE, scale = FALSE)
    object$colsum <- apply(value, 2, sd)
    if (!is.null(object$pCCA) && object$pCCA$rank > 0) {
        comm <- qr.resid(object$pCCA$QR, value)
    }
    if (!is.null(object$CCA) && object$CCA$rank > 0) {
        v <- crossprod(value, object$CCA$u)
        v <- decostand(v, "normalize", MARGIN = 2)
        object$CCA$v <- v
        value <- qr.resid(object$CCA$QR, value)
    }
    if (!is.null(object$CA) && object$CA$rank > 0) {
        v <- crossprod(value, object$CA$u)
        v <- decostand(v, "normalize", MARGIN = 2)
        object$CA$v <- v
    }
    object
}

## capscale may have species scores, but is otherwise similar to dbrda

`sppscores<-.capscale` <-
    function(object, value)
{
    object <- `sppscores<-.dbrda`(object, value)
    object$vdata <- deparse(substitute(value))
    object
}

## metaMDS

`sppscores<-.metaMDS` <-
    function(object, value)
{
    wa <- wascores(object$points, value, expand = TRUE)
    attr(wa, "data") <- deparse(substitute(value))
    object$species <- wa
    object
}

## the main purpose of accessor function is to provide nicer command
## autocompletion and cross-references in help, and of course, to tell
## that it is not implemented (and may never be)

`sppscores` <-
    function(object)
{
    .NotYetImplemented()
}

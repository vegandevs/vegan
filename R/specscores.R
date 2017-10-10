#' Add Species Scores to Ordination Results
#'
#' @param object Ordination object
#' @param comm Community data
#'
#' @rdname specscores
#' @export
`specscores` <-
    function(object, comm)
{
    UseMethod("specscores")
}
#' importFrom stat
#'
#' @rdname specscores
#' @export
`specscores.dbrda` <-
    function(object, comm)
{
    comm <- scale(comm, center = TRUE, scale = FALSE)
    object$colsum <- apply(comm, 2, sd)
    if (!is.null(object$pCCA) && object$pCCA$rank > 0) {
        comm <- qr.resid(object$pCCA$QR, comm)
    }
    if (!is.null(object$CCA) && object$CCA$rank > 0) {
        v <- crossprod(comm, object$CCA$u)
        v <- decostand(v, "normalize", MARGIN = 2)
        object$CCA$v <- v
        comm <- qr.resid(object$CCA$QR, comm)
    }
    if (!is.null(object$CA) && object$CA$rank > 0) {
        v <- crossprod(comm, object$CA$u)
        v <- decostand(v, "normalize", MARGIN = 2)
        object$CA$v <- v
    }
    object
}

## capscale may have species scores, but is otherwise similar to dbrda

`specscores.capscale` <-
    function(object, comm)
{
    if (any(!is.na(object$colsum)))
        warning("function overwrites old species scores")
    specscores.dbrda(object, comm)
}

## metaMDS

`specscores.metaMDS` <-
    function(object, comm, expand = TRUE)
{
    if (any(!is.na(object$species)))
        warning("function overwrites old species scores")
    wa <- wascores(object$points, comm, expand = expand)
    object$species <- wa
    object
}

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
    ## pCCA step looks redundant: see comments in commit d915763d
    if (!is.null(object$pCCA) && object$pCCA$rank > 0) {
        value <- qr.resid(object$pCCA$QR, value)
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

## monoMDS

`sppscores<-.monoMDS` <-
    function(object, value)
{
    wa <- wascores(object$points, value, expand = TRUE)
    attr(wa, "data") <- deparse(substitute(value))
    object$species <- wa
    object
}


## sppscores<-.wcmdscale works with the simple case of row
## weights. Verified to give the same result as RDA with duplicated
## rows when using the number of duplicates as row weights in
## wcmdscale of Eucliden distances (per sqrt(n-1)). However, does not
## give the same results as CA in wcmdscale(vegdist(x, "chisq"),
## w=rowSums(x)/sum(x)).

`sppscores<-.wcmdscale` <-
    function(object, value)
{
    value <- as.matrix(value)
    w <- weights(object)
    spp <- crossprod(value, w * object$points)
    spp <- decostand(spp, "normalize", MARGIN = 2)
    attr(spp, "vdata") <- deparse(substitute(value))
    object$species <- spp
    object
}

## This code would give CA species scores for
##    wcmdscale(vegdist(x, "chisq"), w = rowSums(x)/sum(x), eig=TRUE)
## or CA as weighted PCoA of Chi-square distances

## `sppscores<-.wcmdscale` <-
##     function(object, value)
## {
##     value <- as.matrix(value)
##     value <- decostand(value, "total", MARGIN=2) # for sum(w) = 1
##     spp <- crossprod(value, object$points) # wa: sum(w*x) with sum(w)=1
##     spp <- sweep(spp, 2, object$eig, "/") # for wascores(..., expand=T)
##     object$species <- spp
##     object
## }
## the main purpose of accessor function is to provide nicer command
## autocompletion and cross-references in help, and of course, to tell
## that it is not implemented (and may never be)

`sppscores` <-
    function(object)
{
    .NotYetImplemented()
}

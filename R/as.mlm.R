`as.mlm` <-
function(x) {
    if (is.null(x$CCA))
        stop("'as.mlm' can be used only for constrained ordination")
    UseMethod("as.mlm")
}

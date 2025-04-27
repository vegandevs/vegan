## this thing creates an environment
## the whole point is to create all possible inputs for
## commsim functions only once and reuse them as necessary
## also helps keeping track of updating process for sequential algorithms
## method$mode can be evaluated and use storage mode accordingly
`nullmodel` <-
    function(x, method)
{
    x <- as.matrix(x)
    if (is.null(dim(x)) || length(dim(x)) != 2L)
        stop("'x' must be a matrix-like object")
    ## see PR #742
    if (prod(dim(x)) == 0)
        stop("'x' cannot have NULL rows or columns")
    if (any(is.na(x)))
        stop("'NA' values not allowed")
    if (any(x<0))
        stop("negative values not allowed")
    if (all(x == 0))
        warning("it makes no sense to model all-zero data")
    method <- make.commsim(method)
    if (method$binary)
        x <- ifelse(x > 0, 1L, 0L)
    int <- method$mode == "integer"
    if (int && abs(sum(x) - sum(as.integer(x))) > 10^-6)
        stop("non-integer values not allowed")
    if (int)
        x <- round(x, 0) # round up to closest integer
    storage.mode(x) <- method$mode
    out <- list(
        data=x,
        nrow=as.integer(dim(x)[1L]),
        ncol=as.integer(dim(x)[2L]),
        rowSums=rowSums(x),
        colSums=colSums(x),
        rowFreq=as.integer(rowSums(x > 0)),
        colFreq=as.integer(colSums(x > 0)),
        totalSum=ifelse(int, as.integer(sum(x)), as.double(sum(x))),
        fill=as.integer(sum(x > 0)),
        commsim=method,
        state=if (method$isSeq) x else NULL,
        iter=if (method$isSeq) as.integer(0L) else NULL
        )
#    storage.mode(out$x) <- method$mode
    storage.mode(out$rowSums) <- method$mode
    storage.mode(out$colSums) <- method$mode
    out <- list2env(out, parent=emptyenv())
    class(out) <- c("nullmodel", "environment")
#    class(out) <- "nullmodel"
    out
}

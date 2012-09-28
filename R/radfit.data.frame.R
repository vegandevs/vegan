`radfit.data.frame` <-
    function(x, ...)
{
    ## x *must* have rownames
    rownames(x) <- rownames(x, do.NULL = TRUE)
    out <- apply(x, 1, radfit, ...)
    if (length(out) == 1)
        out <- out[[1]]
    else {
        Call <- match.call()
        class(out) <- "radfit.frame"
    }
    out
}

`radfit.matrix` <-
    function(x, ...)
{
    radfit(as.data.frame(x), ...)
}

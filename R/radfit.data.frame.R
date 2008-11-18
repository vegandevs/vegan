"radfit.data.frame" <-
    function(df, ...)
{
    ## df *must* have rownames
    rownames(df) <- rownames(df, do.NULL = TRUE)
    out <- apply(df, 1, radfit, ...)
    if (length(out) == 1)
        out <- out[[1]]
    else {
        Call <- match.call()
        class(out) <- "radfit.frame"
    }
    out
}

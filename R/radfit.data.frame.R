`radfit.data.frame` <-
    function(x, ...)
{
    ## x *must* have rownames
    rownames(x) <- rownames(x, do.NULL = TRUE)
    ## remove empty rows with no species
    nspec <- specnumber(x)
    if (any(nspec == 0)) {
        warning("removed empty rows with no species")
        x <- x[nspec>0,, drop=FALSE]
    }
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

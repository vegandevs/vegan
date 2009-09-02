### A pair of functions to handle na.action = na.exclude in cca and
### rda (and capscale in the future?). Function ordiNAexclude finds
### the WA scores for NA constraints if possible, and puts these into
### ordination object. Function ordiNApredict pads the result scores
### with NA or scores if available.

`ordiNAexclude` <-
    function(x, excluded)
{
    ## Check that there is a na.action of class "exclude"
    nas <- x$na.action
    if (is.null(nas))
        return(x)
    ## Estimate WA scores for NA cases with newdata of excluded
    ## observations
    if (is.null(x$pCCA)) {
        if (!is.null(x$CCA))
            x$CCA$wa.excluded <- predict(x, newdata = excluded,
                                         type = "wa", model = "CCA")
        if (!is.null(x$CA))
            x$CA$u.excluded <- predict(x, newdata = excluded,
                                       type = "wa", model = "CA")
    }
    x
}

### Put NA or fitted WA among the scores

`ordiNApredict` <-
    function(omit, x)
{
    ## Only do this if omit is of class "exclude"
    if (!inherits(omit, "exclude"))
        return(x)
    x$rowsum <- napredict(omit, x$rowsum) # or zero here?
    if (!is.null(x$CCA)) {
        x$CCA$u <- napredict(omit, x$CCA$u)
        x$CCA$wa <- napredict(omit, x$CCA$wa)
        if (!is.null(x$CCA$wa.excluded))
            x$CCA$wa[omit,] <- x$CCA$wa.excluded
    }
    if (!is.null(x$CA)) {
        x$CA$u <- napredict(omit, x$CA$u)
        if (!is.null(x$CA$u.excluded))
            x$CA$u[omit,] <- x$CA$u.excluded
    }
    x
}

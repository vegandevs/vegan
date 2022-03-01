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
    ## add a 'residuals' item, because step, add1.default and
    ## drop1.default use this to check that number of observations
    ## does not change in sequential fits.
    x$residuals.zombie <- rep(TRUE, max(0, nrow(x$CA$u)))
    ## rowsums for CA (in RDA/PCA rowsum = NA)
    if (!inherits(x, "rda"))
        x$rowsum.excluded <- rowSums(excluded)/x$grand.total
    ## Estimate WA scores for NA cases with newdata of excluded
    ## observations
    if (is.null(x$pCCA) && inherits(nas, "exclude") &&
        !inherits(x, c("dbrda", "capscale"))) {
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
    if (!inherits(x, "rda")) {
        x$rowsum <- napredict(omit, x$rowsum)
        if (inherits(omit, "exclude"))
            x$rowsum[omit] <- x$rowsum.excluded
    }
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

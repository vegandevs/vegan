##' @title Utility for handling user friendly scaling --- None exported
##'
##' @description Convert user-friendly descriptions of scalings to numeric codes used by \code{scores} to date.
##'
##' @param scaling character; which type of scaling is required?
##' @param correlation logical; should correlation-like scores be returned?
##' @param hill logical; should Hill's scaling scores be returned?
`scalingType` <- function(scaling = c("none", "sites", "species", "symmetric"),
                          correlation = FALSE, hill = FALSE)
{
    ## numeric scaling: return as-is
    if (is.numeric(scaling))
        return(scaling)
    ## non-numeric scaling: change to numeric
    tab <- c("none", "sites", "species", "symmetric")
    scaling <- match.arg(scaling)
    scl <- match(scaling, tab) - 1      # -1 as none == scaling 0
    if (scl > 0 && (correlation || hill)) {
        scl <- -scl
    }
    scl                                 # return
}

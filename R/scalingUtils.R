##' @title Utility for handling user friendly scaling --- None exported
##'
##' @description Convert user-friendly descriptions of scalings to numeric codes used by \code{scores} to date.
##'
##' @param scaling character or numeric; which type of scaling is required? Numeric values are returned unaltered
##' @param correlation logical; should correlation-like scores be returned?
##' @param hill logical; should Hill's scaling scores be returned?
`scalingType` <- function(scaling = c("none", "sites", "species", "symmetric"),
                          correlation = FALSE, hill = FALSE) {
    ## Only process scaling further if it is character
    if (is.character(scaling)) {
        ## non-numeric scaling: change to numeric
        tab <- c("none", "sites", "species", "symmetric")
        scaling <- match.arg(scaling)
        scaling <- match(scaling, tab) - 1      # -1 as none == scaling 0
        if (scaling > 0 && (correlation || hill)) {
            scaling <- -scaling
        }
    } else {
        if (!is.numeric(scaling)) {
            stop("'scaling' is neither 'numeric' nor 'character'.")
        }
    }
    scaling                             # return
}

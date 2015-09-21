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
    if (is.numeric(scaling)) {
        return(scaling)                 # numeric; return early
    } else if (is.character(scaling)) {
        ## non-numeric scaling: change to correct numeric code
        scaling <- match.arg(scaling)   # match user choice
        ## Keep `tab` as this is the order of numeric codes
        ## Allows potential to change the default ordering of formal argument 'scaling'
        tab <- c("none", "sites", "species", "symmetric")
        scaling <- match(scaling, tab) - 1      # -1 as none == scaling 0
        if (correlation || hill) {
            scaling <- -scaling
        }
    } else {
        stop("'scaling' is not 'numeric' nor 'character'.")
    }
    scaling                             # return
}

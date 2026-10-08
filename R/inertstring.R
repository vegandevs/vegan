### Get inertia name string for distance-based methods that use
### arguments sqrt.dist and/or add

#' @param inertia Basename of inertia, typically attr(<dist>, "method")
#' @param sqrt.dist passed from function call for taking sqrt(<dist>)
#' @param adjust mean adjustment; in non-adjusted data this is set 1
#' @param ac additive constant returned from addLingoes/Cailliez
#' @param add name of additive adjustment, set to string when evaluating ac
#'
`inertstring` <-
    function(inertia, sqrt.dist, adjust, ac, add)
{
    if (is.null(inertia))
        inertia <- "user-supplied"
    inertia <- paste0(toupper(substr(inertia, 1, 1)),
                     substring(inertia, 2))
    inertia <- paste(inertia, "distance")
    if (!sqrt.dist)
        inertia <- paste("squared", inertia)
    if (adjust != 1) # set to 1, so exact
        inertia <- paste("mean", inertia)
    if (!missing(ac) && ac > sqrt(.Machine$double.eps)) {
        if (is.logical(add) && add)
            add <- "lingoes"
        inertia <- paste(paste0(toupper(substring(add, 1, 1)),
                                substring(add, 2)),
                         "adjusted", inertia)
    }
    inertia
}

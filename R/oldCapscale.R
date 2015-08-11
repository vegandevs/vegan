### Internal function to remove changes in capscale in vegan
### 2.4-0. From vegan 2.4-0 inertia components include negative
### eigenvalues of PCoA whereas earlier ordination and the components
### were only based on real components of PCoA. This function puts
### back the inertia components that ignore imaginary
### dimensions. F-statistics etc will change as a result.  The
### function is provided to maintain compatibility with the old vegan,
### and may be eliminated in the future.
`oldCapscale` <-
    function(object)
{
    ## no imaginary component: nothing need be done
    if (is.null(object$CA$imaginary.rank))
        return(object)
    ## inertia components based only on real dimensions
    object$tot.chi <- object$real.tot.chi
    if (!is.null(object$pCCA)) {
        object$pCCA$tot.chi <- object$pCCA$real.tot.chi
    }
    if (!is.null(object$CCA)) {
        object$CCA$tot.chi <- object$CCA$real.tot.chi
    }
    if (!is.null(object$CA)) {
        object$CA$tot.chi <- object$CA$real.tot.chi
    }
    ## tell what you did
    message("imaginary variation was discarded")
    class(object) <- c("oldcapscale", class(object))
    object
}

# Extract eigenvalues from an object that has them

`eigenvals` <-
    function(x)
{
    UseMethod("eigenvals")
}

`eigenvals.default`<-
    function(x)
{
    ## svd and eigen return unspecified 'list', see if this could be
    ## either of them
    if (is.list(x)) {
        ## eigen
        if (length(x) == 2 && all(names(x) %in% c("values", "vectors")))
            return(x$values)
        ## svd: return squares of singular values
        if (length(x) == 3 && all(names(x) %in% c("d", "u", "v")))
            return(x$d^2)
        else
            stop("No eigenvalues found")
    }
    else
        stop("No eigenvalues found")
}

## squares of sdev 
`eigenvals.prcomp` <-
    function(x)
{
    x$sdev^2
}

## squares of sdev
`eigenvals.princomp` <-
    function(x)
{
    x$sdev^2
}

## concatenate constrained and unconstrained eigenvalues in cca, rda
## and capscale (vegan) -- ignore pCCA component
`eigenvals.cca` <- function(x)
{
   c(x$CCA$eig, x$CA$eig) 
}

## wcmdscale (in vegan)
`eigenvals.wcmdscale` <-
    function(x)
{
    x$eig
}

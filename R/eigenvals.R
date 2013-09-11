# Extract eigenvalues from an object that has them

`eigenvals` <-
    function(x, ...)
{
    UseMethod("eigenvals")
}

`eigenvals.default`<-
    function(x, ...)
{
    ## svd and eigen return unspecified 'list', see if this could be
    ## either of them (like does cmdscale)
    out <- NA
    if (is.list(x)) {
        ## eigen
        if (length(x) == 2 && all(names(x) %in% c("values", "vectors")))
            out <- x$values
        ## svd: return squares of singular values
        else if (length(x) == 3 && all(names(x) %in% c("d", "u", "v")))
            out <- x$d^2
        ## cmdscale() will return all eigenvalues from R  2.12.1
        else if (all(c("points","eig","GOF") %in% names(x)))
            out <- x$eig
    }
    class(out) <- "eigenvals"
    out
}

## squares of sdev
`eigenvals.prcomp` <-
    function(x, ...)
{
    out <- x$sdev^2
    names(out) <- colnames(x$rotation)
    class(out) <- "eigenvals"
    out
}

## squares of sdev
`eigenvals.princomp` <-
    function(x, ...)
{
    out <- x$sdev^2
    class(out) <- "eigenvals"
    out
}

## concatenate constrained and unconstrained eigenvalues in cca, rda
## and capscale (vegan) -- ignore pCCA component
`eigenvals.cca` <- function(x, constrained = FALSE, ...)
{
   if (constrained)
       out <- x$CCA$eig
   else
       out <- c(x$CCA$eig, x$CA$eig)
   if (!is.null(out))
       class(out) <- "eigenvals"
   out
}

## wcmdscale (in vegan)
`eigenvals.wcmdscale` <-
    function(x, ...)
{
    out <- x$eig
    class(out) <- "eigenvals"
    out
}

## pcnm (in vegan)
`eigenvals.pcnm` <-
    function(x, ...)
{
    out <- x$values
    class(out) <- "eigenvals"
    out
}

## betadisper (vegan)
`eigenvals.betadisper` <-  function(x, ...) {
    out <- x$eig
    class(out) <- "eigenvals"
    out
}

## dudi objects of ade4

`eigenvals.dudi` <-
    function(x, ...)
{
    out <- x$eig
    class(out) <- "eigenvals"
    out
}

## labdsv::pco

`eigenvals.pco` <-
    function(x, ...)
{
    out <- x$eig
    class(out) <- "eigenvals"
    out
}

## labdsv::pca

`eigenvals.pca` <-
    function(x, ...)
{
    out <- x$sdev^2
    ## pca() may return only some first eigenvalues
    if ((seig <- sum(out)) < x$totdev) {
        names(out) <- paste("PC", seq_along(out), sep="")
        out <- c(out, "Rest" = x$totdev - seig)
    }
    class(out) <- "eigenvals"
    out
}


`print.eigenvals` <-
    function(x, ...)
{
    print(zapsmall(unclass(x), ...))
    invisible(x)
}

`summary.eigenvals` <-
    function(object, ...)
{
    ## abs(object) is to handle neg eigenvalues of wcmdscale and
    ## capscale
    vars <- object/sum(abs(object))
    importance <- rbind(`Eigenvalue` = object,
                        `Proportion Explained` = round(abs(vars), 5),
                        `Cumulative Proportion`= round(cumsum(abs(vars)), 5))
    out <- list(importance = importance)
    class(out) <- c("summary.eigenvals", "summary.prcomp")
    out
}

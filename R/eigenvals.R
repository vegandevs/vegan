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
    ## honour prcomp(..., rank.=) which only requests rank. eigenvalues
    if (ncol(x$rotation) < length(out)) {
        sumev <- sum(out)
        out <- out[seq_len(ncol(x$rotation))]
        attr(out, "sumev") <- sumev
    }
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
`eigenvals.cca` <- function(x, model = c("all", "unconstrained", "constrained"),
                            constrained = NULL, ...)
{
    out <- if (!is.null(constrained)) {
        ## old behaviour
        message("Argument `constrained` is deprecated; use `model` instead.")
        if (constrained) {
            x$CCA$eig
        } else {
            c(x$CCA$eig, x$CA$eig)
        }
    } else {
        ## new behaviour
        model <- match.arg(model)
        if (identical(model, "all")) {
            c(x$CCA$eig, x$CA$eig)
        } else if (identical(model, "unconstrained")) {
            x$CA$eig
        } else {
            x$CCA$eig
        }
    }
    if (!is.null(out)) {
        class(out) <- "eigenvals"
    }
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

`summary.eigenvals` <- function(object, ...) {
    ## dbRDA can have negative eigenvalues: do not give cumulative
    ## proportions
    if (!is.null(attr(object, "sumev"))) {
        sumev <- attr(object, "sumev")
    } else {
        sumev <- sum(object)
    }
    vars <- object/sumev
    cumvars <- if (all(vars >= 0)) {
        cumsum(vars)
    } else {
        NA
    }
    out <- rbind(`Eigenvalue` = object,
                 `Proportion Explained` = abs(vars),
                 `Cumulative Proportion` = cumvars)
    class(out) <- c("summary.eigenvals", "matrix")
    out
}

## before R svn commit 70391 we used print.summary.prcomp, but now we
## need our own version that is similar to pre-70391 R function

`print.summary.eigenvals` <-
    function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    cat("Importance of components:\n")
    class(x) <- "matrix"
    print(x, digits = digits, ...)
    invisible(x)
}

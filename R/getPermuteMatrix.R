### Interface to the permute package

### input can be (1) a single number giving the number of
### permutations, (2) a how() structure for control parameter in
### permute::shuffleSet, or (3) a permutation matrix which is returned
### as is. In addition, there can be a 'strata' argument which will
### modify case (1). The number of shuffled items must be given in 'N'.

`getPermuteMatrix` <-
    function(perm, N, strata = NULL, model = NULL)
{
    ## Expanded plot assignments can repeat source row indices.
    if (inherits(perm, "plotPermutations")) {
        if (!identical(model, "constraints"))
            stop('plot variable permutations require model = "constraints"')
        if (!isTRUE(ncol(perm) == N))
            stop("plot permutations do not match the fitted observations")
        return(perm)
    }
    if (inherits(perm, "how") && !is.null(getPlots(perm)$variables) &&
        !identical(model, "constraints"))
        stop('plot variable permutations require model = "constraints"')
    if (length(perm) == 1L && is.null(dim(perm)))
        perm <- how(nperm = perm)
    ## apply 'strata', but only if possible: ignore silently other cases
    if (!missing(strata) && !is.null(strata)) {
        if (inherits(perm, "how") && is.null(getBlocks(perm)))
            setBlocks(perm) <- strata
    }
    ## now 'perm' is either a how() or a matrix
    if (inherits(perm, "how"))
        perm <- shuffleSet(N, control = perm)
    else {
        ## The matrix may store integers as doubles.
        if (!is.integer(perm) &&
            (!is.numeric(perm) || is.complex(perm) || anyNA(perm) ||
             any(perm != round(perm))))
            stop("permutation matrix must be strictly integers: use round()")
    }
    ## now 'perm' is a matrix (or always was). If it is a plain
    ## matrix, set minimal attributes for printing. This is a dirty
    ## kluge: should be handled more cleanly.
    if (is.null(attr(perm, "control"))) {
        ## Use valid permute classes instead of a minimal list. Recent
        ## permute versions dispatch getType() on the Within component,
        ## whereas older releases also accept this complete structure.
        control <- how(nperm = nrow(perm))
        control$within$type <- "supplied matrix"
        attr(perm, "control") <- control
    }
    perm
}

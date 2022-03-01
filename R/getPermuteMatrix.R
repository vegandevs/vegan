### Interface to the permute package

### input can be (1) a single number giving the number of
### permutations, (2) a how() structure for control parameter in
### permute::shuffleSet, or (3) a permutation matrix which is returned
### as is. In addition, there can be a 'strata' argument which will
### modify case (1). The number of shuffled items must be given in 'N'.

`getPermuteMatrix` <-
    function(perm, N,  strata = NULL)
{
    ## 'perm' is either a single number, a how() structure or a
    ## permutation matrix
    if (length(perm) == 1) {
        perm <- how(nperm = perm)
    }
    ## apply 'strata', but only if possible: ignore silently other cases
    if (!missing(strata) && !is.null(strata)) {
        if (inherits(perm, "how") && is.null(getBlocks(perm)))
            setBlocks(perm) <- strata
    }
    ## now 'perm' is either a how() or a matrix
    if (inherits(perm, "how"))
        perm <- shuffleSet(N, control = perm)
    else { # matrix: check that it *strictly* integer
        if(!is.integer(perm) && !all(perm == round(perm)))
           stop("permutation matrix must be strictly integers: use round()")
    }
    ## now 'perm' is a matrix (or always was). If it is a plain
    ## matrix, set minimal attributes for printing. This is a dirty
    ## kluge: should be handled more cleanly.
    if (is.null(attr(perm, "control")))
        attr(perm, "control") <-
            structure(list(within=list(type="supplied matrix"),
                           nperm = nrow(perm)), class = "how")
    perm
}

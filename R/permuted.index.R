"permuted.index" <-
    function (n, strata) 
{
    .Deprecated("permute package (shuffle or shuffleSet)")
    if (missing(strata) || is.null(strata)) 
        out <- sample.int(n, n)
    else {
        out <- 1:n
        inds <- names(table(strata))
        for (is in inds) {
            gr <- out[strata == is]
            if (length(gr) > 1) 
                out[gr] <- sample(gr, length(gr))
        }
    }
    out
}


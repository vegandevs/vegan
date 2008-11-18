permute <- function(i, n, control) {
    if(control$complete && !is.null(control$all.perms))
        perm <- control$all.perms[i,]
    else {
        if(control$complete)
            warning("'$all.perms' is NULL, yet '$complete = TRUE'.\nReturning a random permutation.")
        perm <- permuted.index2(n, control)
    }
    return(perm)
}

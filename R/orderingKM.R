`orderingKM` <-
    function(mat)
{

### INPUT :
### mat		(n x k): n the objects and k the descriptors
### 			This matrix must be integers and numeric
###			And must not be binary, it is the partition matrix
###			output by cascadeKM
### OUTPUT :  Ordered matrix

    ## Uses the standard R function cmdscale with matching coefficient
    ## provided by vegdist.c as method=50.

    ##Check up

    if(!is.matrix(mat)) stop("'mat' must be a matrix")
    if(!is.numeric(mat)) stop("'mat' must be numeric")
    if(any(is.na(mat))) stop("'NA' value was found in the matrix")
    if(any(is.infinite(mat))) stop("'Inf' value was found in the matrix")
    nb.desc=ncol(mat)
    nb.obj=nrow(mat)

    d <- .Call(do_vegdist, as.matrix(mat), as.integer(50), PACKAGE = "vegan")
    attributes(d) <- list("Size" = nb.obj,
                          "Labels" = dimnames(mat)[[1]],
                          "Diag" = FALSE,
                          "Upper" = FALSE,
                          "method" = "matching",
                          "class" = "dist")
    scores <- cmdscale(d, k = 1)[,1]
    scores <- order(scores)
    mat[scores,]
}

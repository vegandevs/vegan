`orderingKM` <-
    function(mat)
{

### INPUT :
### mat		(n x k): n the objects and k the descriptors
### 			This matrix must be integers and numeric
###			And must not be binairy, it is the partition matrix
###			output by cascadeKM
### OUTPUT :  Ordered matrix

    ## Uses alternatively the fast USEPOWERALGORITHM provided by Legendre
    ## et al., or the standard R function cmdscale with matching
    ## coefficient provided by vegdist.c as method=50.

    USEPOWERALGORITHM <- FALSE

    ##Check up

    if(!is.matrix(mat)) stop("'mat' must be a matrix")
    if(!is.numeric(mat)) stop("'mat' must be numeric")
    if(any(is.na(mat))) stop("'NA' value was found in the matrix")
    if(any(is.infinite(mat))) stop("'Inf' value was found in the matrix")
    nb.desc=ncol(mat)
    nb.obj=nrow(mat)

    scores<-rep(0.0,nb.obj)
    if (USEPOWERALGORITHM) {
	scores <- as.vector(.Fortran(orderdata, as.integer(mat),
                                 as.integer(nb.obj), as.integer(nb.desc),
                                 sc=as.double(scores), PACKAGE = "vegan")$sc)
    } else {
        d <- .Call(do_vegdist, as.matrix(mat), as.integer(50), PACKAGE = "vegan")
        attributes(d) <- list("Size" = nb.obj,
                              "Labels" = dimnames(mat)[[1]],
                              "Diag" = FALSE,
                              "Upper" = FALSE,
                              "method" = "matching",
                              "class" = "dist")
        scores <- cmdscale(d, k = 1)[,1]
    }
    scores <- order(scores)
    mat[scores,]
}

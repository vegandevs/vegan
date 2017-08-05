`avgdist` <- 
    function(x, sample, distfun = vegdist, meanfun = mean, transf = NULL, iterations = 100, dmethod = "bray", ...)
{
    if (is.na(sample))
        stop("invalid subsampling depth")
    if (is.na(iterations))
        stop("invalid iteration count")
    inputcast <- x
    distfun <- match.fun(distfun)
    # Get the list of iteration matrices
    distlist <- lapply(c(1:iterations), function(i) {
        inputcast <- rrarefy(inputcast, sample = sample)
        # Remove those that did not meet the depth cutoff
        inputcast <- inputcast[c(rowSums(inputcast) %in% sample),]
        if (!is.null(transf)) {
            transmat <- match.fun(transf)
            inputcast <- transmat(inputcast)
        }
        outdist <- distfun(inputcast, method = dmethod, diag = TRUE, upper = TRUE, ...)
        return(as.matrix(outdist))
    })
    # Use the dist list to get the average values
    meanfun <- match.fun(meanfun)
    afunc <- array(
        unlist(as.matrix(distlist)), c(dim(as.matrix(distlist[[1]])), length(distlist)))
    output <- apply(afunc, 1:2, meanfun, ...)
    as.dist(output, upper = TRUE, diag = TRUE)
}

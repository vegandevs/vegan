`avgdist` <- 
    function(x, sample, distfun = vegdist, iterations = 100, dmethod = "bray", ...)
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
        outdist <- distfun(inputcast, method = dmethod, diag = TRUE, upper = TRUE, ...)
        return(outdist)
    })
    # Use the dist list to get the average values
    findist <- Reduce("+", distlist) / length(distlist)
    return(as.dist(findist))
}

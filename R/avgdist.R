`avgdist` <-
    function(x, sample, distfun = vegdist, meanfun = mean,
             transf = NULL, iterations = 100, dmethod = "bray", ...)
{
    if (is.na(sample))
        stop("invalid subsampling depth")
    if (is.na(iterations))
        stop("invalid iteration count")
    distfun <- match.fun(distfun)
    if (!is.null(transf))
        transf <- match.fun(transf)
    # Get the list of iteration matrices
    distlist <- lapply(seq_len(iterations), function(i) {
        inputcast <- rrarefy(x, sample = sample)
        # Remove those that did not meet the depth cutoff
        inputcast <- inputcast[c(rowSums(inputcast) %in% sample),]
        if (!is.null(transf)) {
            inputcast <- transf(inputcast)
        }
        outdist <- distfun(inputcast, method = dmethod,
                           diag = TRUE, upper = TRUE, ...)
        as.matrix(outdist)
    })
    # Use the dist list to get the average values
    meanfun <- match.fun(meanfun)
    afunc <- array(
        unlist(as.matrix(distlist)),
               c(dim(as.matrix(distlist[[1]])), length(distlist)))
    output <- apply(afunc, 1:2, meanfun, ...)
    as.dist(output, upper = TRUE, diag = TRUE)
}

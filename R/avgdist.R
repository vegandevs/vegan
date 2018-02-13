`avgdist` <-
    function(x, sample, distfun = vegdist, meanfun = mean,
             transf = NULL, iterations = 100, dmethod = "bray", ...)
{
    if (is.na(sample))
        stop("invalid subsampling depth")
    if (is.na(iterations))
        stop("invalid iteration count")
    inputcast <- x
    distfun <- match.fun(distfun)
    if (!is.null(transf))
        transf <- match.fun(transf)
    # Get the list of iteration matrices
    distlist <- lapply(seq_len(iterations), function(i) {
        # Suppress warnings because it will otherwise return many warnings about
        # subsampling depth not being met, which we deal with below by returning
        # samples that do not meet the threshold.
        inputcast <- suppressWarnings(rrarefy(inputcast, sample = sample))
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
    # Save row names from distlist
    # Take from first element since should all be the same
    rnames <- row.names(distlist[[1]])
    afunc <- array(
        unlist(as.matrix(distlist)),
               c(dim(as.matrix(distlist[[1]])), length(distlist)))
    output <- apply(afunc, 1:2, meanfun, ...)
    # Set the names on the matrix
    colnames(output) <- rownames(output) <- rnames
    # Print any samples that were removed, if they were removed
    if(nrow(x) != nrow(output)) {
        dropsamples <- setdiff(row.names(inputcast), row.names(output))
        warning(
            gettextf(
                "The following sampling units were removed because they were below sampling depth: %s",
                paste(dropsamples, collapse = ", ")))
    }
    output <- as.dist(output, diag = TRUE, upper = TRUE)
    attr(output, "call") <- match.call()
    attr(output, "method") <- "avgdist"
    output
}

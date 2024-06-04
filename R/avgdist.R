`avgdist` <-
    function(x, sample, distfun = vegdist, meanfun = mean,
             transf = NULL, iterations = 100, dmethod = "bray",
             diag = TRUE, upper = TRUE, ...)
{
    if (missing(sample)) {
        stop("Subsampling depth must be supplied via argument 'sample'")
    } else {
        if (!(is.numeric(sample) && sample > 0L)) {
            stop("Invalid subsampling depth; 'sample' must be positive & numeric")
        }
    }
    if (!is.numeric(iterations)) {
        stop("Invalid iteration count; must be numeric")
    }
    inputcast <- x
    distfun <- match.fun(distfun)
    if (!is.null(transf)) {
        transf <- match.fun(transf)
    }
    ## warn here if data do not look observed counts with singletons
    minobs <- min(x[x > 0])
    if (minobs > 1)
        warning(gettextf("most observed count data have counts 1, but smallest count is %d", minobs))
    # Get the list of iteration matrices
    distlist <- lapply(seq_len(iterations), function(i) {
        # Suppress warnings because it will otherwise return many warnings about
        # subsampling depth not being met, which we deal with below by returning
        # samples that do not meet the threshold.
        inputcast <- suppressWarnings(rrarefy(inputcast, sample = sample))
        # Remove those that did not meet the depth cutoff
        inputcast <- inputcast[c(rowSums(inputcast) >= sample), ]
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
    dropsamples <- setdiff(row.names(inputcast), row.names(output))
    if (length(dropsamples) > 0L) {
        warning(gettextf(
            "The following sampling units were removed because they were below sampling depth: %s",
                         paste(dropsamples, collapse = ", ")))
    }
    output <- as.dist(output, diag = diag, upper = upper)
    attr(output, "call") <- match.call()
    attr(output, "method") <- paste("avgdist", dmethod)
    output
}

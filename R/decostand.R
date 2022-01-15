`decostand` <-
    function (x, method, MARGIN, range.global, logbase = 2, na.rm = FALSE, ...)
{
    wasDataFrame <- is.data.frame(x)
    x <- as.matrix(x)
    METHODS <- c("total", "max", "frequency", "normalize", "range", "rank",
                 "rrank", "standardize", "pa", "chi.square", "hellinger",
                 "log", "clr", "rclr", "alr", "ilr")
    method <- match.arg(method, METHODS)
    if (any(x < 0, na.rm = na.rm)) {
        k <- min(x, na.rm = na.rm)
        if (method %in% c("total", "frequency", "pa", "chi.square", "rank",
                          "rrank", "clr", "rclr", "alr", "ilr")) {
            warning("input data contains negative entries: result may be non-sense\n")
        }
    }
    else k <- .Machine$double.eps
    switch(method, total = {
        if (missing(MARGIN))
            MARGIN <- 1
        tmp <- pmax(k, apply(x, MARGIN, sum, na.rm = na.rm))
        x <- sweep(x, MARGIN, tmp, "/")
    }, max = {
        if (missing(MARGIN))
            MARGIN <- 2
        tmp <- pmax(k, apply(x, MARGIN, max, na.rm = na.rm))
        x <- sweep(x, MARGIN, tmp, "/")
    }, frequency = {
        if (missing(MARGIN))
            MARGIN <- 2
        tmp <- pmax(k, apply(x, MARGIN, sum, na.rm = na.rm))
        fre <- apply(x > 0, MARGIN, sum, na.rm = na.rm)
        tmp <- fre/tmp
        x <- sweep(x, MARGIN, tmp, "*")
    }, normalize = {
        if (missing(MARGIN))
            MARGIN <- 1
        tmp <- apply(x^2, MARGIN, sum, na.rm = na.rm)
        tmp <- pmax(.Machine$double.eps, sqrt(tmp))
        x <- sweep(x, MARGIN, tmp, "/")
    }, range = {
        if (missing(MARGIN))
            MARGIN <- 2
        if (missing(range.global))
            xtmp <- x
        else {
            if (dim(range.global)[MARGIN] != dim(x)[MARGIN])
                stop("range matrix does not match data matrix")
            xtmp <- as.matrix(range.global)
        }
        tmp <- apply(xtmp, MARGIN, min, na.rm = na.rm)
        ran <- apply(xtmp, MARGIN, max, na.rm = na.rm)
        ran <- ran - tmp
        ran <- pmax(k, ran, na.rm = na.rm)
        x <- sweep(x, MARGIN, tmp, "-")
        x <- sweep(x, MARGIN, ran, "/")
    }, rank = {
        if (missing(MARGIN)) MARGIN <- 1
        x[x==0] <- NA
        x <- apply(x, MARGIN, rank, na.last = "keep")
        if (MARGIN == 1) # gives transposed x
            x <- t(x)
        x[is.na(x)] <- 0
    }, rrank = {
        if (missing(MARGIN)) MARGIN <- 1
        x <- decostand(x, "rank", MARGIN = MARGIN)
        x <- sweep(x, MARGIN, specnumber(x, MARGIN = MARGIN), "/")
    }, standardize = {
        if (!missing(MARGIN) && MARGIN == 1)
            x <- t(scale(t(x)))
        else x <- scale(x)
    }, pa = {
        x <- ifelse(x > 0, 1, 0)
    }, chi.square = {
        if (!missing(MARGIN) && MARGIN == 2)
            x <- t(x)
        x <- sqrt(sum(x, na.rm = na.rm)) * x/outer(pmax(k, rowSums(x,
                         na.rm = na.rm)), sqrt(colSums(x, na.rm = na.rm)))
    }, hellinger = {
        x <- sqrt(decostand(x, "total", MARGIN = MARGIN, na.rm = na.rm))
    }, log = {### Marti Anderson logs, after Etienne Laliberte
        if (!isTRUE(all.equal(as.integer(x), as.vector(x)))) {
            x <- x / min(x[x > 0], na.rm = TRUE)
            warning("non-integer data: divided by smallest positive value",
                    call. = FALSE)
        }
        x[x > 0 & !is.na(x)] <- log(x[x > 0 & !is.na(x)], base = logbase) + 1

    }, alr = {
        if (missing(MARGIN))
	    MARGIN <- 1
        if (MARGIN == 1) 
          x <- .calc_alr(x, ...)
	else x <- t(.calc_alr(t(x), ...))
    }, ilr = {
        if (missing(MARGIN))
	    MARGIN <- 1
        if (MARGIN == 1) 
          x <- .calc_ilr(x, ...)
	else x <- t(.calc_ilr(t(x), ...))
    }, clr = {
        if (missing(MARGIN))
	    MARGIN <- 1
        if (MARGIN == 1) 
          x <- .calc_clr(x, ...)
	else x <- t(.calc_clr(t(x), ...))
    }, rclr = {
        if (missing(MARGIN))
	    MARGIN <- 1
        if (MARGIN == 1) 
          x <- .calc_rclr(x, ...)
	else x <- t(.calc_rclr(t(x), ...))
    })
    if (any(is.nan(x)))
        warning("result contains NaN, perhaps due to impossible mathematical 
                 operation\n")
    if (wasDataFrame)
        x <- as.data.frame(x)
    attr(x, "decostand") <- method
    x
}


# Modified from the original version in mia R package
.calc_clr <- function(x, pseudocount=0, na.rm=TRUE){
    # Add pseudocount
    x <- x + pseudocount
    # If there are negative values, gives an error.
    if (any(x <= 0, na.rm = TRUE)) {
        stop("Abundance table contains zero or negative values and ",
             "clr-transformation is being applied without (suitable) ",
             "pseudocount. \n")
    }
    # In every sample, calculates the log of individual entries.
    # After that calculates
    # the sample-specific mean value and subtracts every entries'
    # value with that.
    #clog <- t(log(x))
    #t(clog - rowMeans(clog))

    clog <- log(x)
    clog - rowMeans(clog)
    
}

# Modified from the original version in mia R package
.calc_rclr <- function(x, na.rm=TRUE){
    # If there are negative values, gives an error.
    if (any(x < 0, na.rm = na.rm)) {
        stop("Abundance table contains negative values. The 
              rclr transformation assumes non-negative values.\n")
    }
   # Log transform
   clog <- log(x)
   # zeros are converted into infinite values in clr
   # They are converted to NAs for now
   clog[is.infinite(clog)] <- NA
   # Calculates means for every sample, does not take NAs into account
   mean_clog <- rowMeans(clog, na.rm = TRUE)
   # Calculates exponential values from means, i.e., geometric means
   geometric_means_of_samples <- exp(mean_clog)
   # Divides all values by their sample-wide geometric means
   # Then does logarithmic transform and transposes the table back to its original
   # form
   xx <- log(x/geometric_means_of_samples)
   # If there were zeros, there are infinite values after logarithmic transform.
   # They are converted to zero.
   xx[is.infinite(xx)] <- 0
   xx
}


.calc_alr <- function (x, reference = 1, na.rm=TRUE, pseudocount=0) {
    # Add pseudocount
    x <- x + pseudocount
    # If there is negative values, gives an error.
    if (any(x < 0, na.rm = na.rm)) {
        stop("Abundance table contains negative values and ",
             "alr-transformation is being applied without (suitable) ",
             "pseudocount. \n")
    }    
    if (reference > nrow(x)) 
        stop("The reference should be a feature name, or index between 1 to", ncol(x))
    clog <- log(x)
    clog[, -reference]-clog[, reference]
}



.calc_ilr <- function (x, pseudocount=0) {
    # Add pseudocount
    x <- x + pseudocount
    # If there is negative values, gives an error.
    if (any(x < 0, na.rm = TRUE)) {
        stop("Abundance table contains negative values and ",
             "alr-transformation is being applied without (suitable) ",
             "pseudocount. \n")
    }    

    # For a more efficient implementation ideas, check the packages
    # compositions and philr for ilrBase
    x.ilr <- matrix(NA, nrow(x), ncol(x)-1)
    rownames(x.ilr) <- rownames(x)    
    for (i in seq_len(nrow(x))) {
        for (j in seq_len(ncol(x.ilr))) {
            x.ilr[i, j] <- -sqrt(j/(j + 1)) * log(((prod(x[i, seq_len(j)]))^(1/j))/x[i, j + 1])
        }
    }
    
    x.ilr

}



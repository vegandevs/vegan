"decostand" <-
    function (x, method, MARGIN, range.global, logbase = 2, na.rm = FALSE, ...) 
{
    wasDataFrame <- is.data.frame(x)
    x <- as.matrix(x)
    METHODS <- c("total", "max", "frequency", "normalize", "range", 
                 "standardize", "pa", "chi.square", "hellinger", "log")
    method <- match.arg(method, METHODS)
    if (any(x < 0, na.rm = na.rm)) {
        k <- min(x, na.rm = na.rm)
        if (method %in% c("total", "frequency", "pa", "chi.square")) {
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
        tmp <- pmax(k, sqrt(tmp))
        x <- sweep(x, MARGIN, tmp, "/")
    }, range = {
        if (missing(MARGIN)) 
            MARGIN <- 2
        if (missing(range.global)) 
            xtmp <- x
        else {
            if (dim(range.global)[MARGIN] != dim(x)[MARGIN]) 
                stop("range matrix doesn't match data matrix")
            xtmp <- as.matrix(range.global)
        }
        tmp <- apply(xtmp, MARGIN, min, na.rm = na.rm)
        ran <- apply(xtmp, MARGIN, max, na.rm = na.rm)
        ran <- ran - tmp
        ran <- pmax(k, ran, na.rm = na.rm)
        x <- sweep(x, MARGIN, tmp, "-")
        x <- sweep(x, MARGIN, ran, "/")
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
    })
    if (any(is.nan(x))) 
        warning("result contains NaN, perhaps due to impossible mathematical operation\n")
    if (wasDataFrame)
        x <- as.data.frame(x)
    attr(x, "decostand") <- method
    x
}

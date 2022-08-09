`decostand` <-
    function (x, method, MARGIN, range.global, logbase = 2, na.rm = FALSE, ...)
{
    wasDataFrame <- is.data.frame(x)
    x <- as.matrix(x)
    METHODS <- c("total", "max", "frequency", "normalize", "range", "rank",
                 "rrank", "standardize", "pa", "chi.square", "hellinger",
                 "log", "clr", "rclr", "alr")
    method <- match.arg(method, METHODS)
    if (any(x < 0, na.rm = na.rm)) {
        k <- min(x, na.rm = na.rm)
        if (method %in% c("total", "frequency", "pa", "chi.square", "rank",
                          "rrank", "clr", "rclr")) {
            warning("input data contains negative entries: result may be non-sense")
        }
    }
    else k <- .Machine$double.eps
    attr <- NULL
    switch(method, total = {
        if (missing(MARGIN))
            MARGIN <- 1
        tmp <- pmax(k, apply(x, MARGIN, sum, na.rm = na.rm))
        x <- sweep(x, MARGIN, tmp, "/")
        attr <- list("total" = tmp, "margin" = MARGIN)
    }, max = {
        if (missing(MARGIN))
            MARGIN <- 2
        tmp <- pmax(k, apply(x, MARGIN, max, na.rm = na.rm))
        x <- sweep(x, MARGIN, tmp, "/")
        attr <- list("max" = tmp, "margin" = MARGIN)
    }, frequency = {
        if (missing(MARGIN))
            MARGIN <- 2
        tmp <- pmax(k, apply(x, MARGIN, sum, na.rm = na.rm))
        fre <- apply(x > 0, MARGIN, sum, na.rm = na.rm)
        tmp <- fre/tmp
        x <- sweep(x, MARGIN, tmp, "*")
        attr <- list("scale" = tmp, "margin" = MARGIN)
    }, normalize = {
        if (missing(MARGIN))
            MARGIN <- 1
        tmp <- apply(x^2, MARGIN, sum, na.rm = na.rm)
        tmp <- pmax(.Machine$double.eps, sqrt(tmp))
        x <- sweep(x, MARGIN, tmp, "/")
        attr <- list("norm" = tmp, "margin" = MARGIN)
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
        attr <- list("min" = tmp, "range" = ran, "margin" = MARGIN)
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
        else {
            x <- scale(x)
            MARGIN <- 2
        }
        attr <- list("center" = attr(x, "scaled:center"),
                     "scale" = attr(x, "scaled:scale"),
                     "margin" = MARGIN)
    }, pa = {
        x <- ifelse(x > 0, 1, 0)
    }, chi.square = {
        if (missing(MARGIN))
            MARGIN <- 1
        ## MARGIN 2 transposes the result!
        if (MARGIN == 2)
            x <- t(x)
        rs <- pmax(k, rowSums(x, na.rm = na.rm))
        cs <- pmax(k, colSums(x, na.rm = na.rm))
        tot <- sum(x, na.rm = na.rm)
        x <- sqrt(tot) * x/outer(rs, sqrt(cs))
        attr <- list("tot" = tot, "rsum" = rs, "csum" = cs, margin = MARGIN)
    }, hellinger = {
        x <- sqrt(decostand(x, "total", MARGIN = MARGIN, na.rm = na.rm))
        attr <- attr(x, "parameters")
    }, log = {### Marti Anderson logs, after Etienne Laliberte
        if (!isTRUE(all.equal(as.integer(x), as.vector(x)))) {
            minpos <- min(x[x > 0], na.rm = TRUE)
            x <- x / minpos
            warning("non-integer data: divided by smallest positive value",
                    call. = FALSE)
        } else {
            minpos <- 1
        }
        x[x > 0 & !is.na(x)] <- log(x[x > 0 & !is.na(x)], base = logbase) + 1
        attr <- list("logbase" = logbase, minpos = minpos)
    }, alr = {
        if (missing(MARGIN))
	    MARGIN <- 1
        if (MARGIN == 1)
            x <- t(.calc_alr(t(x), ...))
	else x <- .calc_alr(x, ...)
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
    attr(x, "parameters") <- attr
    attr(x, "decostand") <- method
    x
}


# Modified from the original version in mia R package
.calc_clr <- function(x, pseudocount=0, na.rm=TRUE){
    # Add pseudocount
    x <- x + pseudocount
    # Error with negative values
    if (any(x <= 0, na.rm = TRUE)) {
        stop("Abundance table contains zero or negative values and ",
             "clr-transformation is being applied without (suitable) ",
             "pseudocount. \n")
    }
    # In every sample, calculate the log of individual entries.
    # Then calculate
    # the sample-specific mean value and subtract every entries'
    # value with that.
    clog <- log(x)
    clog - rowMeans(clog)
}

# Modified from the original version in mia R package
.calc_rclr <- function(x, na.rm=TRUE){
    # Error with negative values
    if (any(x < 0, na.rm = na.rm)) {
        stop("Matrix has negative values but the
              rclr transformation assumes non-negative values.\n")
    }
   # Log transform
   clog <- log(x)
   # Convert zeros to NAs in rclr
   clog[is.infinite(clog)] <- NA
   # Calculate mean for every sample, ignoring the NAs
   mean_clog <- rowMeans(clog, na.rm = TRUE)
   # Calculate geometric means per sample
   geom_mean <- exp(mean_clog)
   # Divide all values by their sample-wide geometric means
   # Log and transpose back to original shape
   xx <- log(x/geom_mean)
   # If there were zeros, there are infinite values after logarithmic transform.
   # Convert those to zero.
   xx[is.infinite(xx)] <- 0
   xx
}


.calc_alr <- function (x, reference = 1, na.rm = TRUE, pseudocount = 0)
{
    # Add pseudocount
    x <- x + pseudocount
    # If there is negative values, gives an error.
    if (any(x < 0, na.rm = na.rm)) {
        stop("cannot be used with negative data: add pseudocount > ",
             abs(max(x[x < 0], na.rm = na.rm)), call. = FALSE)
    }
    ## name must be changed to numeric index for [-reference,] to work
    if (is.character(reference)) {
        reference <- which(reference == colnames(x))
        if (!length(reference)) # found it?
            stop("'reference' name was not found in data", call. = FALSE)
    }
    if (reference > ncol(x) || reference < 1)
        stop("'reference' should be a name or index 1 to ",
             ncol(x), call. = FALSE)
    clog <- log(x)
    clog[, -reference] - clog[, reference]
}

`decobackstand` <-
    function(x, zap = TRUE)
{
    method <- attr(x, "decostand")
    if (is.null(method))
        stop("function can be used only with 'decostand' standardized data")
    para <- attr(x, "parameters")
    x <- switch(method,
                "total" = sweep(x, para$margin, para$total, "*"),
                "max" = sweep(x, para$margin, para$max, "*"),
                "frequency" = sweep(x, para$margin, para$scale, "/"),
                "normalize" = sweep(x, para$margin, para$norm, "*"),
                "range" = { x <- sweep(x, para$margin, para$range, "*")
                            sweep(x, para$margin, para$min, "+")},
                "standardize" = {x <- sweep(x, para$margin, para$scale, "*")
                                 sweep(x, para$margin, para$center, "+") },
                "hellinger" = sweep(x^2, para$margin, para$total, "*"),
                "chi.square" = { rc <- outer(para$rsum, sqrt(para$csum))
                                 x <- x * rc /sqrt(para$tot)
                                 if (para$margin == 1) x else t(x) },
                "log" = { x[x > 0 & !is.na(x)] <-
                              para$logbase^(x[x > 0 & !is.na(x)] - 1)
                              x * para$minpos},
                "wisconsin" = { x <- sweep(x, 1, para$total, "*")
                                sweep(x, 2, para$max, "*") },
                stop("no back-transformation available for method ",
                     sQuote(method))
                )
    if (zap)
        x[abs(x) < sqrt(.Machine$double.eps)] <- 0
    x
}

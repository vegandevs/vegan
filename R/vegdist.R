`vegdist` <-
    function (x, method = "bray", binary = FALSE, diag = FALSE, upper = FALSE,
              na.rm = FALSE, ...)
{
    ZAP <- 1e-15
    if (!is.na(pmatch(method, "euclidian")))
        method <- "euclidean"
    ## the order of METHODS below *MUST* match the #define'd numbers
    ## in vegdist.c
    METHODS <- c("manhattan", "euclidean", "canberra", "bray", # 4
                 "kulczynski", "gower", "morisita", "horn", #8
                 "mountford", "jaccard", "raup", "binomial", "chao", #13
                 "altGower", "cao", "mahalanobis", "clark", "chisq", "chord", #19
		 "hellinger", "aitchison", "robust.aitchison") # 22
    method <- pmatch(method, METHODS)
    inm <- METHODS[method]
    if (is.na(method))
        stop("invalid distance method")
    if (method == -1)
        stop("ambiguous distance method")
    ## most tests are faster for matrix than for data frame, and we
    ## need matrix in .Call() anyway
    x <- as.matrix(x)
    if (!na.rm && anyNA(x))
        stop("missing values are not allowed with argument 'na.rm = FALSE'")
    ## all vegdist indices need numeric data (Gower included).
    if (!(is.numeric(x) || is.logical(x)))
        stop("input data must be numeric")
    if (!method %in% c(1,2,6,16,18) && any(rowSums(x, na.rm = TRUE) == 0))
        warning("you have empty rows: their dissimilarities may be
                 meaningless in method ",
                 dQuote(inm))
    ## 1 manhattan, 2 euclidean, 3 canberra, 6 gower, 16 mahalanobis, 19 chord
    if (!method %in% c(1,2,3,6,16,19,20) && any(x < 0, na.rm = TRUE))
        warning("results may be meaningless because data have negative entries
                 in method ",
                 dQuote(inm))
    if (method %in% c(11,18) && any(colSums(x) == 0, na.rm = TRUE))
        warning("data have empty species which influence the results in
                 method ",
                dQuote(inm))
    if (method == 6) # gower, but no altGower
        x <- decostand(x, "range", 2, na.rm = TRUE, ...)
    if (method == 16) # mahalanobis
        x <- veganMahatrans(scale(x, scale = FALSE), na.rm = na.rm)
    if (method == 18) # chisq
        x <- decostand(x, "chi.square", na.rm = na.rm)
    if (method == 21)  # aitchison
        x <- decostand(x, "clr", ...)  # dots to pass possible pseudocount
    if (method == 22)  # robust.aitchison
        x <- decostand(x, "rclr", na.rm = na.rm, ...) # No pseudocount for rclr
    if (binary)
        x <- decostand(x, "pa")
    N <- nrow(x)
    if (method %in% c(7, 13, 15) && !isTRUE(all.equal(x, round(x))))
        warning("results may be meaningless with non-integer data in method ",
                dQuote(inm))
    d <- .Call(do_vegdist, x, as.integer(method))
    d[d < ZAP] <- 0
    if (any(is.na(d)))
        warning("missing values in results")
    ## add attribute maxdist: the maximum value of the distance function
    attr(d, "maxdist") <-
        if(method %in% c(3,4,5,7,8,10,11,13,17)) # index in 0..1
            1
        else if (method %in% c(19,20)) # chord, hellinger
            sqrt(2)
        else if (method == 9) # Mountford
            log(2)
        else # no fixed upper limit
            NA
    attr(d, "Size") <- N
    attr(d, "Labels") <- dimnames(x)[[1]]
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- paste(if (binary)
                               "binary ", METHODS[method], sep = "")
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    d
}

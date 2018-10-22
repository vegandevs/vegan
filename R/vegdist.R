`vegdist` <-
    function (x, method = "bray", binary = FALSE, diag = FALSE, upper = FALSE,
              na.rm = FALSE, ...)
{
    ZAP <- 1e-15
    if (!is.na(pmatch(method, "euclidian")))
        method <- "euclidean"
    ## the order of METHODS below *MUST* match the #define'd numbers
    ## in vegdist.c
    METHODS <- c("manhattan", "euclidean", "canberra", "bray",
                 "kulczynski", "gower", "morisita", "horn",
                 "mountford", "jaccard", "raup", "binomial", "chao",
                 "altGower", "cao", "mahalanobis", "clark")
    method <- pmatch(method, METHODS)
    inm <- METHODS[method]
    if (is.na(method))
        stop("invalid distance method")
    if (method == -1)
        stop("ambiguous distance method")
    ## most tests are faster for matrix than for data frame, and we
    ## need matrix in .Call() anyway
    x <- as.matrix(x)
    ## all vegdist indices need numeric data (Gower included).
    if (!(is.numeric(x) || is.logical(x)))
        stop("input data must be numeric")
    if (!method %in% c(1,2,6,16) && any(rowSums(x, na.rm = TRUE) == 0))
        warning("you have empty rows: their dissimilarities may be meaningless in method ",
                dQuote(inm))
    ## 1 manhattan, 2 euclidean, 3 canberra, 6 gower, 16 mahalanobis
    if (!method %in% c(1,2,3,6,16) && any(x < 0, na.rm = TRUE))
        warning("results may be meaningless because data have negative entries in method ",
                dQuote(inm))
    if (method == 11 && any(colSums(x) == 0))
        warning("data have empty species which influence the results in method ",
                dQuote(inm))
    if (method == 6) # gower, but no altGower
        x <- decostand(x, "range", 2, na.rm = TRUE, ...)
    if (method == 16) # mahalanobis
        x <- veganMahatrans(scale(x, scale = FALSE))
    if (binary)
        x <- decostand(x, "pa")
    N <- nrow(x)
    if (method %in% c(7, 13, 15) && !identical(all.equal(x, round(x)), TRUE))
        warning("results may be meaningless with non-integer data in method ",
                dQuote(inm))
    d <- .Call(do_vegdist, x, as.integer(method))
    if (method == 10)
        d <- 2 * d/(1 + d)
    d[d < ZAP] <- 0
    if (any(is.na(d)))
        warning("missing values in results")
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

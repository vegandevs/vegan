`vegdist` <-
    function (x, method = "bray", binary = FALSE, diag = FALSE, upper = FALSE, 
              na.rm = FALSE, ...) 
{
    ZAP <- 1e-15
    if (!is.na(pmatch(method, "euclidian"))) 
        method <- "euclidean"
    METHODS <- c("manhattan", "euclidean", "canberra", "bray", 
                 "kulczynski", "gower", "morisita", "horn", "mountford", 
                 "jaccard", "raup", "binomial", "chao", "altGower")
    method <- pmatch(method, METHODS)
    inm <- METHODS[method]
    if (is.na(method)) 
        stop("invalid distance method")
    if (method == -1) 
        stop("ambiguous distance method")
    if (method > 2 && any(rowSums(x, na.rm = TRUE) == 0)) 
        warning("you have empty rows: their dissimilarities may be meaningless in method ", inm,"\n")
    if (method > 2 && any(x < 0, na.rm = TRUE)) 
        warning("results may be meaningless because data have negative entries in method ", inm,"\n")
    if (method == 11 && any(colSums(x) == 0)) 
        warning("data have empty species which influence the results im method ", inm, "\n")
    if (method == 6) # gower, but no altGower
        x <- decostand(x, "range", 2, na.rm = TRUE, ...)
    if (binary) 
        x <- decostand(x, "pa")
    N <- nrow(x <- as.matrix(x))
    if (method %in% c(7, 13) && !identical(all.equal(as.integer(x), 
                                                     as.vector(x)), TRUE)) 
        warning("results may be meaningless with non-integer data in method ", inm, "\n")
    d <- .C("veg_distance", x = as.double(x), nr = N, nc = ncol(x), 
            d = double(N * (N - 1)/2), diag = as.integer(FALSE), 
            method = as.integer(method), NAOK = na.rm, PACKAGE = "vegan")$d
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

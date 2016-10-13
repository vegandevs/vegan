`renyi` <-
    function (x, scales = c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64,
                 Inf), hill = FALSE)
{
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    if (p == 1) {
        x <- t(x)
        n <- nrow(x)
        p <- ncol(x)
    }
    ## do not make total=1 if not needed (diversity() does anyway,
    ## species richness does not need)
    if (!all(scales %in% c(0,1)))
        x <- sweep(x, 1, rowSums(x), "/")
    m <- length(scales)
    result <- array(0, dim = c(n, m))
    dimnames(result) <- list(sites = rownames(x), scale = scales)
    for (a in 1:m) {
        if (scales[a] != 0 && scales[a] != 1 && scales[a] !=  Inf) {
            result[, a] <- log(rowSums(x^scales[a]))/(1 - scales[a])
        }
        else {
            if (scales[a] == 0) {
                result[, a] <- log(rowSums(x > 0))
            }
            else if (scales[a] == Inf) {
                result[, a] <- -log(apply(x, 1, max))
            }
            else {
                result[, a] <- diversity(x)
            }
        }
    }
    if (hill)
        result <- exp(result)
    if (any(dim(result) == 1))
        result <- drop(result)
    else
        result <- as.data.frame(result)
    class(result) <- c("renyi", class(result))
    result
}

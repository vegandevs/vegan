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
    x <- decostand(x, "total", 1)
    m <- length(scales)
    result <- array(0, dim = c(n, m))
    dimnames(result) <- list(sites = rownames(x), scale = scales)
    for (a in 1:m) {
        if (scales[a] != 0 && scales[a] != 1 && scales[a] != 
            Inf) {
            result[, a] <- log(apply(x^scales[a], 1, sum))/(1 - 
                                                            scales[a])
        }
        else {
            if (scales[a] == 0) {
                result[, a] <- log(apply(x > 0, 1, sum))
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
    result <- as.data.frame(result)
    if (any(dim(result) == 1)) 
        result <- unlist(result, use.names = TRUE)
    class(result) <- c("renyi", class(result))
    result
}

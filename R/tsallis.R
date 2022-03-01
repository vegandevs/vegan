`tsallis` <-
    function (x, scales = seq(0, 2, 0.2), norm=FALSE, hill=FALSE)
{
    if (norm && hill)
        stop("'norm = TRUE' and 'hill = TRUE' should not be used at the same time")
    x <- as.matrix(x)
    if (!is.numeric(x))
        stop("input data must be numeric")
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
        if (scales[a] != 1 && scales[a] != 0) {
                result[, a] <- (1-(rowSums(x^scales[a])))/(scales[a] - 1)
        }
        else {
            if (scales[a] == 1) result[, a] <- diversity(x, "shannon")
            if (scales[a] == 0) result[, a] <- rowSums(x > 0) - 1
        }
        if (norm) {
            ST <- rowSums(x > 0)
            if (scales[a] == 1) result[, a] <- result[, a] / log(ST)
            else result[, a] <- result[, a] / ((ST^(1-scales[a]) - 1) / (1 - scales[a]))
        }
        if (hill) {
            result[, a] <- if (scales[a] == 1) {
                exp(result[, a])
            } else {
                (1 - (scales[a] - 1) * result[, a])^(1/(1-scales[a]))
            }
        }
    }
    result <- as.data.frame(result)
    if (any(dim(result) == 1))
        result <- unlist(result, use.names = TRUE)
    class(result) <- c("tsallis", "renyi", class(result))
    result
}

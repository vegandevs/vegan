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
    ## scale rows to unit total
    x <- sweep(x, 1, rowSums(x), "/")
    m <- length(scales)
    result <- array(0, dim = c(n, m))
    dimnames(result) <- list(sites = rownames(x), scale = scales)
    for (a in 1:m) {
        result[,a] <-
            switch(as.character(scales[a]),
                   "0" = log(rowSums(x > 0)),
                   "1" = -rowSums(x * log(x), na.rm = TRUE),
                   "2" = -log(rowSums(x^2)),
                   "Inf" =  -log(apply(x, 1, max)),
                   log(rowSums(x^scales[a]))/(1 - scales[a]))
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

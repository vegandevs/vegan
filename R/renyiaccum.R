`renyiaccum` <-
function(x, scales=c(0, 0.5, 1, 2, 4, Inf), permutations = 100,
         raw = FALSE, collector = FALSE, subset, ...)
{
    if (!missing(subset))
        x <- subset(x, subset)
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    if (p==1) {
        x <- t(x)
        n <- nrow(x)
        p <- ncol(x)
    }
    pmat <- getPermuteMatrix(permutations, n)
    permutations <- nrow(pmat)
    m <- length(scales)
    result <- array(dim=c(n,m,permutations))
    dimnames(result) <- list(pooled.sites=c(1:n), scale=scales,
                             permutation=c(1:permutations))
    for (k in 1:permutations) {
        result[,,k] <- as.matrix(renyi((apply(x[pmat[k,],],2,cumsum)),
                                       scales=scales, ...))
    }
    if (raw)
        collector <- FALSE
    if (collector)
        ref <- as.matrix(renyi(apply(x, 2, cumsum), scales = scales, ...))
    if (raw) {
        if (m==1) {
            result <- result[,1,]
        }
    }else{
        tmp <- array(dim=c(n,m,6 + as.numeric(collector)))
        for (i in 1:n) {
            for (j in 1:m) {
                tmp[i,j,1] <- mean(result[i,j,1:permutations]) 
                tmp[i,j,2] <- sd(result[i,j,1:permutations])
                tmp[i,j,3] <- min(result[i,j,1:permutations])
                tmp[i,j,4] <- max(result[i,j,1:permutations])
                tmp[i,j,5] <- quantile(result[i,j,1:permutations],0.025)
                tmp[i,j,6] <- quantile(result[i,j,1:permutations],0.975)
                if (collector)
                    tmp[i,j,7] <- ref[i,j]
            }
        }
        result <- tmp
        dimnames(result) <- list(pooled.sites=c(1:n),
                                  scale=scales,
                                  c("mean", "stdev", "min", "max", "Qnt 0.025", "Qnt 0.975", if (collector) "Collector"))
    }
    attr(result, "control") <- attr(pmat, "control")
    class(result) <- c("renyiaccum", class(result))
    result
}


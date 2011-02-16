## Contribution diversity
## Lu, H.P., H.H. Wagner and X.Y. Chen (2007). 
## A contribution diversity approach to evaluate species diversity. 
## Basic and Applied Ecology 8: 1 -12.
`contribdiv` <-
    function(comm, index = c("richness", "simpson"), relative = FALSE,
             scaled = TRUE, drop.zero = FALSE)
{

    index <- match.arg(index)

    x <- comm[rowSums(comm) > 0, colSums(comm) > 0]
    n <- nrow(x)
    S <- ncol(x)

    if (index == "richness") {
        n.i <- colSums(x > 0)
        S.k <- rowSums(x > 0)
        alpha <- S.k / n
        beta <- apply(x, 1, function(z) sum((n - n.i[z > 0]) / (n * n.i[z > 0])))
        denom <- 1
    } else {
        P.ik <- decostand(x, "total")
        P.i <- apply(P.ik, 2, function(z) sum(z) / n)
        P.i2 <- matrix(P.i, n, S, byrow=TRUE)
        alpha <- diversity(x, "simpson")
        beta <- rowSums(P.ik * (P.ik - P.i2))
        denom <- n
    }
    gamma <- alpha + beta
    D <- sum(beta) / sum(gamma)
    if (relative) {
        denom <- if (scaled)
            {denom * sum(gamma)} else 1
        alpha <- (alpha - mean(alpha)) / denom
        beta <- (beta - mean(beta)) / denom
        gamma <- (gamma - mean(gamma)) / denom
    }
    rval <- data.frame(alpha = alpha, beta = beta, gamma = gamma)
    if (!drop.zero && nrow(comm) != n) {
        nas <- rep(NA, nrow(comm))
        rval2 <- data.frame(alpha = nas, beta = nas, gamma = nas)
        rval2[rowSums(comm) > 0, ] <- rval
        rval <- rval2
    }
    attr(rval, "diff.coef") <- D
    attr(rval, "index") <- index
    attr(rval, "relative") <- relative
    attr(rval, "scaled") <- scaled
    class(rval) <- c("contribdiv", "data.frame")
    rval
}

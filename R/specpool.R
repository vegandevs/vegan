"specpool" <-
    function (x, pool) 
{
    x <- as.matrix(x)
    if (missing(pool)) 
        pool <- rep("All", nrow(x))
    out <- seq(1:nrow(x))
    groups <- table(pool)
    inds <- names(groups)
    S <- var.chao <- chao <- var.jack1 <- jack.1 <- jack.2 <- var.boot <- bootS <-
        rep(NA, length(inds))
    names(S) <- names(var.chao) <- names(chao) <- names(var.jack1) <-
        names(jack.1) <- names(jack.2) <- names(var.boot) <- names(bootS) <- inds
    for (is in inds) {
        a1 <- a2 <- NA
        gr <- out[pool == is]
        n <- length(gr)
        X <- x[gr, , drop = FALSE]
        freq <- colSums(X > 0)
        p <- freq[freq > 0]/n
        S[is] <- sum(freq > 0)
        if (S[is] == 0) 
            next
        if (n >= 1) 
            a1 <- sum(freq == 1)
        if (n >= 2) 
            a2 <- sum(freq == 2)
        else 0
        chao[is] <- S[is] + if(!is.na(a2) && a2 > 0)
            a1 * a1/2/a2
        else 0
        jack.1[is] <- S[is] + a1 * (n - 1)/n
        jack.2[is] <- S[is] + a1 * (2 * n - 3)/n - a2 * (n - 
                                                         2)^2/n/(n - 1)
        bootS[is] <- S[is] + sum((1 - p)^n)
        aa <- if (!is.na(a2) && a2 > 0) 
            a1/a2
        else 0
        var.chao[is] <- a2 * (0.5 + (1 + aa/4) * aa) * aa * aa
        if (!is.na(a1) && a1 > 0) {
            jf <- table(rowSums(X[, freq == 1, drop = FALSE] > 
                                0))
            var.jack1[is] <- (sum(as.numeric(names(jf))^2 * jf) - 
                              a1/n) * (n - 1)/n
        }
        pn <- (1 - p)^n
        X <- X[, freq > 0, drop = FALSE]
        Zp <- (crossprod(X == 0)/n)^n - outer(pn, pn, "*")
        var.boot[is] <- sum(pn * (1 - pn)) + 2 * sum(Zp[lower.tri(Zp)])
    }
    out <- list(Species = S, Chao = chao,  Chao.SE = sqrt(var.chao), 
                Jack.1 = jack.1, Jack1.SE = sqrt(var.jack1), Jack.2 = jack.2, 
                Boot = bootS, Boot.SE = sqrt(var.boot), n = as.vector(groups))
    out <- as.data.frame(out)
    attr(out, "pool") <- pool
    out
}

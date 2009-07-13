`dispindmorisita` <-
function(x, unique.rm=FALSE, crit=0.05)
{
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    Imor <- numeric(p)
    Imor <- apply(x, 2, function(y) n * ((sum(y^2) - sum(y)) / (sum(y)^2 - sum(y))))
    Smor <- Imor
    chicr <- qchisq(c(0+crit/2, 1-crit/2), n-1, lower.tail=FALSE)
    Muni <- apply(x, 2, function(y) (chicr[2] - n + sum(y)) / (sum(y) - 1))
    Mclu <- apply(x, 2, function(y) (chicr[1] - n + sum(y)) / (sum(y) - 1))
    for (i in 1:p) {
        if (apply(x, 2, sum)[i] > 1) {
            if (Imor[i] >= Mclu[i] && Mclu[i] > 1)
                Smor[i] <- 0.5 + 0.5 * ((Imor[i] - Mclu[i]) / (n - Mclu[i]))
            if (Mclu[i] > Imor[i] && Imor[i] >=1)
                Smor[i] <- 0.5 * ((Imor[i] - 1) / (Mclu[i] - 1))
            if (1 > Imor[i] && Imor[i] > Muni[i])
                Smor[i] <- -0.5 * ((Imor[i] - 1) / (Muni[i] - 1))
            if (1 > Muni[i] && Muni[i] > Imor[i])
                Smor[i] <- -0.5 + 0.5 * ((Imor[i] - Muni[i]) / Muni[i])
        }
    }
    out <- data.frame(imor=Imor,mclu=Mclu,muni=Muni,imst=Smor)
    usp <- which(apply(x > 0, 2, sum) == 1)
    if (unique.rm && length(usp) != 0)
        out <- out[-usp,]
    return(out)
}


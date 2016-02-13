"downweight" <-
function (veg, fraction = 5) 
{
    Const1 <- 1e-10
    if (fraction < 1) 
        fraction <- 1/fraction
    veg <- as.matrix(veg)
    yeig1 <- apply(veg, 2, sum)
    y2 <- apply(veg^2, 2, sum) + Const1
    y2 <- yeig1^2/y2
    amax <- max(y2)/fraction
    v <- rep(1, ncol(veg))
    downers <- y2 < amax
    v[downers] <- (y2/amax)[downers]
    veg <- sweep(veg, 2, v, "*")
    attr(veg, "v") <- v
    attr(veg, "fraction") <- fraction
    veg
}

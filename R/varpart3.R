"varpart3" <-
    function (Y, X1, X2, X3) 
{
    Y <- as.matrix(Y)
    X1 <- as.matrix(X1)
    X2 <- as.matrix(X2)
    X3 <- as.matrix(X3)
    n <- nrow(Y)
    n1 <- nrow(X1)
    n2 <- nrow(X2)
    n3 <- nrow(X3)
    p <- ncol(Y)
    mm1 <- ncol(X1)
    mm2 <- ncol(X2)
    mm3 <- ncol(X3)
    if (n1 != n) 
        stop("Y and X1 do not have the same number of rows")
    if (n2 != n) 
        stop("Y and X2 do not have the same number of rows")
    if (n3 != n) 
        stop("Y and X3 do not have the same number of rows")
    Y <- scale(Y, center = TRUE, scale = FALSE)
    X1 <- scale(X1, center = TRUE, scale = TRUE)
    X2 <- scale(X2, center = TRUE, scale = TRUE)
    X3 <- scale(X3, center = TRUE, scale = TRUE)
    SS.Y <- sum(Y * Y)
    dummy <- simpleRDA2(Y, X1, SS.Y, mm1)
    adfg.ua <- dummy$Rsquare
    m1 <- dummy$m
    if (m1 != mm1) 
        warning("collinearity detected in X1: mm = ", mm1, ", m = ", 
                m1, call. = FALSE)
    dummy <- simpleRDA2(Y, X2, SS.Y, mm2)
    bdeg.ua <- dummy$Rsquare
    m2 <- dummy$m
    if (m2 != mm2) 
        warning("collinearity detected in X2: mm = ", mm2, ", m = ", 
                m2, call. = FALSE)
    dummy <- simpleRDA2(Y, X3, SS.Y, mm3)
    cefg.ua <- dummy$Rsquare
    m3 <- dummy$m
    if (m3 != mm3) 
        warning("collinearity detected in X3: mm =", mm3, ", m =", 
                m3, call. = FALSE)
    mm4 = mm1 + mm2
    dummy <- simpleRDA2(Y, cbind(X1, X2), SS.Y, mm4)
    abdefg.ua <- dummy$Rsquare
    m4 <- dummy$m
    if (m4 != mm4) 
        warning("collinearity detected in cbind(X1,X2): mm = ", 
                mm4, ", m = ", m4, call. = FALSE)
    mm5 = mm1 + mm3
    dummy <- simpleRDA2(Y, cbind(X1, X3), SS.Y, mm5)
    acdefg.ua <- dummy$Rsquare
    m5 <- dummy$m
    if (m5 != mm5) 
        warning("collinearity detected in cbind(X1,X3): mm = ", 
                mm5, ", m = ", m5, call. = FALSE)
    mm6 = mm2 + mm3
    dummy <- simpleRDA2(Y, cbind(X2, X3), SS.Y, mm6)
    bcdefg.ua <- dummy$Rsquare
    m6 <- dummy$m
    if (m6 != mm6) 
        warning("collinearity detected in cbind(X2,X3): mm = ", 
                mm6, ", m = ", m6, call. = FALSE)
    mm7 = mm1 + mm2 + mm3
    dummy <- simpleRDA2(Y, cbind(X1, X2, X3), SS.Y, mm7)
    abcdefg.ua <- dummy$Rsquare
    m7 <- dummy$m
    if (m7 != mm7) 
        warning("collinearity detected in cbind(X1,X2,X3): mm = ", 
                mm7, ", m = ", m7, call. = FALSE)
    bigwarning <- NULL
    if ((m1 + m2) > m4) 
        bigwarning <- c(bigwarning, c("X1, X2"))
    if ((m1 + m3) > m5) 
        bigwarning <- c(bigwarning, c("X1, X3"))
    if ((m2 + m3) > m6) 
        bigwarning <- c(bigwarning, c("X2, X3"))
    if ((m1 + m2 + m3) > m7) 
        bigwarning <- c(bigwarning, c("X1, X2, X3"))
    adfg <- RsquareAdj(adfg.ua, n, m1)
    bdeg <- RsquareAdj(bdeg.ua, n, m2)
    cefg <- RsquareAdj(cefg.ua, n, m3)
    abdefg <- RsquareAdj(abdefg.ua, n, m4)
    acdefg <- RsquareAdj(acdefg.ua, n, m5)
    bcdefg <- RsquareAdj(bcdefg.ua, n, m6)
    abcdefg <- RsquareAdj(abcdefg.ua, n, m7)
    Df <- c(m1, m2, m3, m4, m5, m6, m7)
    fract <- data.frame(Df = Df,
                        R.square = c(adfg.ua, bdeg.ua, cefg.ua, abdefg.ua, acdefg.ua,
                        bcdefg.ua, abcdefg.ua), 
                        Adj.R.square = c(adfg, bdeg, cefg, abdefg, acdefg, bcdefg, 
                        abcdefg), Testable = rep(TRUE, 7) & Df)
    rownames(fract) <- c("[a+d+f+g] = X1", "[b+d+e+g] = X2", 
                         "[c+e+f+g] = X3", "[a+b+d+e+f+g] = X1+X2", "[a+c+d+e+f+g] = X1+X3", 
                         "[b+c+d+e+f+g] = X2+X3", "[a+b+c+d+e+f+g] = All")
    a <- abcdefg - bcdefg
    b <- abcdefg - acdefg
    c <- abcdefg - abdefg
    d <- acdefg - cefg - a
    e <- abdefg - adfg - b
    f <- bcdefg - bdeg - c
    g <- adfg - a - d - f
    ma <- m7 - m6
    mb <- m7 - m5
    mc <- m7 - m4
    mad <- m5 - m3
    maf <- m4 - m2
    mbd <- m6 - m3
    mbe <- m4 - m1
    mce <- m5 - m1
    mcf <- m6 - m2
    Df <-  c(ma, mb, mc, rep(0, 4), NA)
    indfract <- data.frame(Df = Df, 
                           R.square = rep(NA, 8),
                           Adj.R.square = c(a, b, c, d, e, f, g, 1 - abcdefg),
                           Testable = c(rep(TRUE, 3), rep(FALSE, 5)) & Df)
    rownames(indfract) <- c("[a] = X1 | X2+X3", "[b] = X2 | X1+X3", 
                            "[c] = X3 | X1+X2", "[d]", "[e]", "[f]", "[g]", "[h] = Residuals")
    Df <- c(mad, maf, mbd, mbe, mce, mcf)
    contr1 <- data.frame(Df = Df,
                         R.square = rep(NA, 6),
                         Adj.R.square = c(a + d, a + f, b + d, b + e, c + e, c + f),
                         Testable = rep(TRUE, 6) & Df)
    rownames(contr1) <- c("[a+d] = X1 | X3", "[a+f] = X1 | X2", 
                          "[b+d] = X2 | X3", "[b+e] = X2 | X1", "[c+e] = X3 | X1", 
                          "[c+f] = X3 | X2")
    out <- list(fract = fract, indfract = indfract, contr1 = contr1, 
                SS.Y = SS.Y, nsets = 3, bigwarning = bigwarning, n = n1)
    class(out) <- "varpart234"
    out
}

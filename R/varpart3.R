`varpart3` <-
    function (Y, X1, X2, X3, chisquare, permat)
{
    collinwarn <- function(case, mm, m)
        warning(gettextf("collinearity detected in %s: mm = %d, m = %d",
                         case, mm, m), call. = FALSE)
    if (inherits(Y, "dist")) {
        Y <- GowerDblcen(as.matrix(Y^2), na.rm = FALSE)
        Y <- -Y/2
        SS.Y <- sum(diag(Y))
        simpleRDA2 <- match.fun(simpleDBRDA)
    } else {
        Y <- as.matrix(Y)
        if (chisquare) {
            SS.Y <- sum(initCA(Y)^2)
            simpleRDA2 <- match.fun(simpleCCA)
        } else {
            Y <- scale(Y, center = TRUE, scale = FALSE)
            SS.Y <- sum(Y * Y)
        }
    }
    X1 <- as.matrix(X1)
    X2 <- as.matrix(X2)
    X3 <- as.matrix(X3)
    n <- nrow(Y)
    n1 <- nrow(X1)
    n2 <- nrow(X2)
    n3 <- nrow(X3)
    mm1 <- ncol(X1)
    mm2 <- ncol(X2)
    mm3 <- ncol(X3)
    if (n1 != n)
        stop("Y and X1 do not have the same number of rows")
    if (n2 != n)
        stop("Y and X2 do not have the same number of rows")
    if (n3 != n)
        stop("Y and X3 do not have the same number of rows")
    X1 <- scale(X1, center = TRUE, scale = FALSE)
    X2 <- scale(X2, center = TRUE, scale = FALSE)
    X3 <- scale(X3, center = TRUE, scale = FALSE)
    dummy <- simpleRDA2(Y, X1, SS.Y, permat)
    adfg.ua <- dummy$Rsquare
    adfg <- dummy$RsquareAdj
    m1 <- dummy$m
    if (m1 != mm1)
        collinwarn("X1", mm1, m1)
    dummy <- simpleRDA2(Y, X2, SS.Y, permat)
    bdeg.ua <- dummy$Rsquare
    bdeg <- dummy$RsquareAdj
    m2 <- dummy$m
    if (m2 != mm2)
        collinwarn("X2", mm2, m2)
    dummy <- simpleRDA2(Y, X3, SS.Y, permat)
    cefg.ua <- dummy$Rsquare
    cefg <- dummy$RsquareAdj
    m3 <- dummy$m
    if (m3 != mm3)
        collinwarn("X3", mm3, m3)
    mm4 = mm1 + mm2
    dummy <- simpleRDA2(Y, cbind(X1, X2), SS.Y, permat)
    abdefg.ua <- dummy$Rsquare
    abdefg <- dummy$RsquareAdj
    m4 <- dummy$m
    if (m4 != mm4)
        collinwarn("cbind(X1,X2)", mm4, m4)
    mm5 = mm1 + mm3
    dummy <- simpleRDA2(Y, cbind(X1, X3), SS.Y, permat)
    acdefg.ua <- dummy$Rsquare
    acdefg <- dummy$RsquareAdj
    m5 <- dummy$m
    if (m5 != mm5)
        collinwarn("cbind(X1,X3)", mm5, m5)
    mm6 = mm2 + mm3
    dummy <- simpleRDA2(Y, cbind(X2, X3), SS.Y, permat)
    bcdefg.ua <- dummy$Rsquare
    bcdefg <- dummy$RsquareAdj
    m6 <- dummy$m
    if (m6 != mm6)
        collinwarn("cbind(X2,X3)", mm6, m6)
    mm7 = mm1 + mm2 + mm3
    dummy <- simpleRDA2(Y, cbind(X1, X2, X3), SS.Y, permat)
    abcdefg.ua <- dummy$Rsquare
    abcdefg <- dummy$RsquareAdj
    m7 <- dummy$m
    if (m7 != mm7)
        collinwarn("cbind(X1,X2,X3)", mm7, m7)
    bigwarning <- NULL
    if ((m1 + m2) > m4)
        bigwarning <- c(bigwarning, c("X1, X2"))
    if ((m1 + m3) > m5)
        bigwarning <- c(bigwarning, c("X1, X3"))
    if ((m2 + m3) > m6)
        bigwarning <- c(bigwarning, c("X2, X3"))
    if ((m1 + m2 + m3) > m7)
        bigwarning <- c(bigwarning, c("X1, X2, X3"))
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

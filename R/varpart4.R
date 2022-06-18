`varpart4` <-
    function (Y, X1, X2, X3, X4, chisquare, permat)
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
            SS.Y <- sum(initCA(Y^2))
            simpleRDA2 <- match.fun(simpleCCA)
        } else {
            Y <- scale(Y, center = TRUE, scale = FALSE)
            SS.Y <- sum(Y * Y)
        }
    }
    X1 <- as.matrix(X1)
    X2 <- as.matrix(X2)
    X3 <- as.matrix(X3)
    X4 <- as.matrix(X4)
    n <- nrow(Y)
    n1 <- nrow(X1)
    n2 <- nrow(X2)
    n3 <- nrow(X3)
    n4 <- nrow(X4)
    mm1 <- ncol(X1)
    mm2 <- ncol(X2)
    mm3 <- ncol(X3)
    mm4 <- ncol(X4)
    if (n1 != n)
        stop("Y and X1 do not have the same number of rows")
    if (n2 != n)
        stop("Y and X2 do not have the same number of rows")
    if (n3 != n)
        stop("Y and X3 do not have the same number of rows")
    if (n4 != n)
        stop("Y and X4 do not have the same number of rows")
    X1 <- scale(X1, center = TRUE, scale = FALSE)
    X2 <- scale(X2, center = TRUE, scale = FALSE)
    X3 <- scale(X3, center = TRUE, scale = FALSE)
    X4 <- scale(X4, center = TRUE, scale = FALSE)
    dummy <- simpleRDA2(Y, X1, SS.Y, permat)
    aeghklno.ua <- dummy$Rsquare
    aeghklno <- dummy$RsquareAdj
    m1 <- dummy$m
    if (m1 != mm1)
        collinwarn("X1", mm1, m1)
    dummy <- simpleRDA2(Y, X2, SS.Y, permat)
    befiklmo.ua <- dummy$Rsquare
    befiklmo <- dummy$RsquareAdj
    m2 <- dummy$m
    if (m2 != mm2)
        collinwarn("X2", mm2, m2)
    dummy <- simpleRDA2(Y, X3, SS.Y, permat)
    cfgjlmno.ua <- dummy$Rsquare
    cfgjlmno <- dummy$RsquareAdj
    m3 <- dummy$m
    if (m3 != mm3)
        collinwarn("X3", mm3, m3)
    dummy <- simpleRDA2(Y, X4, SS.Y, permat)
    dhijkmno.ua <- dummy$Rsquare
    dhijkmno <- dummy$RsquareAdj
    m4 <- dummy$m
    if (m4 != mm4)
        collinwarn("X4", mm4, m4)
    mm5 = mm1 + mm2
    dummy <- simpleRDA2(Y, cbind(X1, X2), SS.Y, permat)
    abefghiklmno.ua <- dummy$Rsquare
    abefghiklmno <- dummy$RsquareAdj
    m5 <- dummy$m
    if (m5 != mm5)
        collinwarn("cbind(X1,X2)", mm5, m5)
    mm6 = mm1 + mm3
    dummy <- simpleRDA2(Y, cbind(X1, X3), SS.Y, permat)
    acefghjklmno.ua <- dummy$Rsquare
    acefghjklmno <- dummy$RsquareAdj
    m6 <- dummy$m
    if (m6 != mm6)
        collinwarn("cbind(X1,X3", mm6, m6)
    mm7 = mm1 + mm4
    dummy <- simpleRDA2(Y, cbind(X1, X4), SS.Y, permat)
    adeghijklmno.ua <- dummy$Rsquare
    adeghijklmno <- dummy$RsquareAdj
    m7 <- dummy$m
    if (m7 != mm7)
        collinwarn("cbind(X1,X4)", mm7, m7)
    mm8 = mm2 + mm3
    dummy <- simpleRDA2(Y, cbind(X2, X3), SS.Y, permat)
    bcefgijklmno.ua <- dummy$Rsquare
    bcefgijklmno <- dummy$RsquareAdj
    m8 <- dummy$m
    if (m8 != mm8)
        collinwarn("cbind(X2,X3)", mm8, m8)
    mm9 = mm2 + mm4
    dummy <- simpleRDA2(Y, cbind(X2, X4), SS.Y, permat)
    bdefhijklmno.ua <- dummy$Rsquare
    bdefhijklmno <- dummy$RsquareAdj
    m9 <- dummy$m
    if (m9 != mm9)
        collinwarn("cbind(X2,X4)", mm9, m9)
    mm10 = mm3 + mm4
    dummy <- simpleRDA2(Y, cbind(X3, X4), SS.Y, permat)
    cdfghijklmno.ua <- dummy$Rsquare
    cdfghijklmno <- dummy$RsquareAdj
    m10 <- dummy$m
    if (m10 != mm10)
        collinwarn("cbind(X3,X4)", mm10, m10)
    mm11 = mm1 + mm2 + mm3
    dummy <- simpleRDA2(Y, cbind(X1, X2, X3), SS.Y, permat)
    abcefghijklmno.ua <- dummy$Rsquare
    abcefghijklmno <- dummy$RsquareAdj
    m11 <- dummy$m
    if (m11 != mm11)
        collinwarn("cbind(X1,X2,X3)", mm11, m11)
    mm12 = mm1 + mm2 + mm4
    dummy <- simpleRDA2(Y, cbind(X1, X2, X4), SS.Y, permat)
    abdefghijklmno.ua <- dummy$Rsquare
    abdefghijklmno <- dummy$RsquareAdj
    m12 <- dummy$m
    if (m12 != mm12)
        collinwarn("c(X1,X2,X4)", mm12, m12)
    mm13 = mm1 + mm3 + mm4
    dummy <- simpleRDA2(Y, cbind(X1, X3, X4), SS.Y, permat)
    acdefghijklmno.ua <- dummy$Rsquare
    acdefghijklmno <- dummy$RsquareAdj
    m13 <- dummy$m
    if (m13 != mm13)
        collinwarn("cbind(X1,X3,X4)", mm13, m13)
    mm14 = mm2 + mm3 + mm4
    dummy <- simpleRDA2(Y, cbind(X2, X3, X4), SS.Y, permat)
    bcdefghijklmno.ua <- dummy$Rsquare
    bcdefghijklmno <- dummy$RsquareAdj
    m14 <- dummy$m
    if (m14 != mm14)
        collinwarn("cbind(X2,X3,X4)", mm14, m14)
    mm15 = mm1 + mm2 + mm3 + mm4
    dummy <- simpleRDA2(Y, cbind(X1, X2, X3, X4), SS.Y, permat)
    abcdefghijklmno.ua <- dummy$Rsquare
    abcdefghijklmno <- dummy$RsquareAdj
    m15 <- dummy$m
    if (m15 != mm15)
        collinwarn("cbind(X1,X2,X3,X4)", mm15, m15)
    bigwarning <- NULL
    if ((m1 + m2) > m5)
        bigwarning <- c(bigwarning, c("X1, X2"))
    if ((m1 + m3) > m6)
        bigwarning <- c(bigwarning, c("X1, X3"))
    if ((m1 + m4) > m7)
        bigwarning <- c(bigwarning, c("X1, X4"))
    if ((m2 + m3) > m8)
        bigwarning <- c(bigwarning, c("X2, X3"))
    if ((m2 + m4) > m9)
        bigwarning <- c(bigwarning, c("X2, X4"))
    if ((m3 + m4) > m10)
        bigwarning <- c(bigwarning, c("X3, X4"))
    if ((m1 + m2 + m3) > m11)
        bigwarning <- c(bigwarning, c("X1, X2, X3"))
    if ((m1 + m2 + m4) > m12)
        bigwarning <- c(bigwarning, c("X1, X2, X4"))
    if ((m1 + m3 + m4) > m13)
        bigwarning <- c(bigwarning, c("X1, X3, X4"))
    if ((m2 + m3 + m4) > m14)
        bigwarning <- c(bigwarning, c("X2, X3, X4"))
    if ((m1 + m2 + m3 + m4) > m15)
        bigwarning <- c(bigwarning, c("X1, X2, X3, X4"))
    Df <-  c(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15)
    fract <- data.frame(Df = Df,
                        R.square = c(aeghklno.ua, befiklmo.ua, cfgjlmno.ua,
                        dhijkmno.ua, abefghiklmno.ua, acefghjklmno.ua,
                        adeghijklmno.ua, bcefgijklmno.ua, bdefhijklmno.ua,
                        cdfghijklmno.ua, abcefghijklmno.ua, abdefghijklmno.ua,
                        acdefghijklmno.ua, bcdefghijklmno.ua, abcdefghijklmno.ua),
                        Adj.R.square = c(aeghklno, befiklmo, cfgjlmno, dhijkmno,
                        abefghiklmno, acefghjklmno, adeghijklmno, bcefgijklmno,
                        bdefhijklmno, cdfghijklmno, abcefghijklmno, abdefghijklmno,
                        acdefghijklmno, bcdefghijklmno, abcdefghijklmno),
                        Testable = rep(TRUE, 15) & Df)
    rownames(fract) <- c("[aeghklno] = X1", "[befiklmo] = X2",
                         "[cfgjlmno] = X3", "[dhijkmno] = X4", "[abefghiklmno] = X1+X2",
                         "[acefghjklmno] = X1+X3", "[adeghijklmno] = X1+X4", "[bcefgijklmno] = X2+X3",
                         "[bdefhijklmno] = X2+X4", "[cdfghijklmno] = X3+X4", "[abcefghijklmno] = X1+X2+X3",
                         "[abdefghijklmno] = X1+X2+X4", "[acdefghijklmno] = X1+X3+X4",
                         "[bcdefghijklmno] = X2+X3+X4", "[abcdefghijklmno] = All")
    ae = acdefghijklmno - cdfghijklmno
    ag = abdefghijklmno - bdefhijklmno
    ah = abcefghijklmno - bcefgijklmno
    be = bcdefghijklmno - cdfghijklmno
    bf = abdefghijklmno - adeghijklmno
    bi = abcefghijklmno - acefghjklmno
    cf = acdefghijklmno - adeghijklmno
    cg = bcdefghijklmno - bdefhijklmno
    cj = abcefghijklmno - abefghiklmno
    dh = bcdefghijklmno - bcefgijklmno
    di = acdefghijklmno - acefghjklmno
    dj = abdefghijklmno - abefghiklmno
    Df <- c(m13-m10, m12-m9, m11-m8, m14-m10, m12-m7, m11-m6, m13-m7, m14-m9,
            m11-m5, m14-m8, m13-m6, m12-m5)
    contr2 <- data.frame(Df = Df,
                         R.square = rep(NA, 12),
                         Adj.R.square = c(ae, ag, ah,
                         be, bf, bi, cf, cg, cj, dh, di, dj),
                         Testable = rep(TRUE, 12) & Df)
    rownames(contr2) <- c("[ae] = X1 | X3+X4", "[ag] = X1 | X2+X4",
                          "[ah] = X1 | X2+X3", "[be] = X2 | X3+X4", "[bf] = X2 | X1+X4",
                          "[bi] = X2 | X1+X3", "[cf] = X3 | X1+X4", "[cg] = X3 | X2+X4",
                          "[cj] = X3 | X1+X2", "[dh] = X4 | X2+X3", "[di] = X4 | X1+X3",
                          "[dj] = X4 | X1+X2")
    aghn = abefghiklmno - befiklmo
    aehk = acefghjklmno - cfgjlmno
    aegl = adeghijklmno - dhijkmno
    bfim = abefghiklmno - aeghklno
    beik = bcefgijklmno - cfgjlmno
    befl = bdefhijklmno - dhijkmno
    cfjm = acefghjklmno - aeghklno
    cgjn = bcefgijklmno - befiklmo
    cfgl = cdfghijklmno - dhijkmno
    dijm = adeghijklmno - aeghklno
    dhjn = bdefhijklmno - befiklmo
    dhik = cdfghijklmno - cfgjlmno
    Df <- c(m5-m2, m6-m3, m7-m4, m5-m1, m8-m3, m9-m4, m6-m1, m8-m2, m10-m4,
            m7-m1, m9-m2, m10-m3)
    contr1 <- data.frame(Df = Df,
                         R.square = rep(NA, 12),
                         Adj.R.square = c(aghn, aehk,
                         aegl, bfim, beik, befl, cfjm, cgjn, cfgl, dijm, dhjn,
                                                 dhik),
                         Testable = rep(TRUE, 12) & Df)
    rownames(contr1) <- c("[aghn] = X1 | X2", "[aehk] = X1 | X3",
                          "[aegl] = X1 | X4", "[bfim] = X2 | X1", "[beik] = X2 | X3",
                          "[befl] = X2 | X4", "[cfjm] = X3 | X1", "[cgjn] = X3 | X2",
                          "[cfgl] = X3 | X4", "[dijm] = X4 | X1 ", "[dhjn] = X4 | X2",
                          "[dhik] = X4 | X3")
    a <- abcdefghijklmno - bcdefghijklmno
    b <- abcdefghijklmno - acdefghijklmno
    c <- abcdefghijklmno - abdefghijklmno
    d <- abcdefghijklmno - abcefghijklmno
    e <- ae - a
    f <- bf - b
    g <- ag - a
    h <- ah - a
    i <- bi - b
    j <- cj - c
    k <- aehk - ae - h
    l <- aegl - ae - g
    m <- bfim - bf - i
    n <- aghn - ag - h
    o <- aeghklno - aehk - g - l - n
    indfract <- data.frame(Df = c(m15-m14, m15-m13, m15-m12, m15-m11, rep(0, 12)),
                           R.square = rep(NA, 16), Adj.R.square = c(a, b, c, d,
                                                   e, f, g, h, i, j, k, l, m, n, o, 1 - abcdefghijklmno),
                           Testable = c(rep(TRUE, 4), rep(FALSE, 12)))
    rownames(indfract) <- c("[a] = X1 | X2+X3+X4", "[b] = X2 | X1+X3+X4",
                            "[c] = X3 | X1+X2+X4", "[d] = X4 | X1+X2+X3", "[e]",
                            "[f]", "[g]", "[h]", "[i]", "[j]", "[k]", "[l]", "[m]",
                            "[n]", "[o]", "[p] = Residuals")
    out <- list(fract = fract, indfract = indfract, contr1 = contr1,
                contr2 = contr2, SS.Y = SS.Y, nsets = 4, bigwarning = bigwarning,
                n = n1)
    class(out) <- "varpart234"
    out
}

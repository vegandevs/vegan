"varpart2" <-
    function (Y, X1, X2) 
{
    Y <- as.matrix(Y)
    X1 <- as.matrix(X1)
    X2 <- as.matrix(X2)
    n <- nrow(Y)
    n1 <- nrow(X1)
    n2 <- nrow(X2)
    p <- ncol(Y)
    mm1 <- ncol(X1)
    mm2 <- ncol(X2)
    if (n1 != n) 
        stop("Y and X1 do not have the same number of rows")
    if (n2 != n) 
        stop("Y and X2 do not have the same number of rows")
    Y <- scale(Y, center = TRUE, scale = FALSE)
    X1 <- scale(X1, center = TRUE, scale = FALSE)
    X2 <- scale(X2, center = TRUE, scale = FALSE)
    SS.Y <- sum(Y * Y)
    dummy <- simpleRDA2(Y, X1, SS.Y, mm1)
    ab.ua <- dummy$Rsquare
    m1 <- dummy$m
    if (m1 != mm1) 
        warning("collinearity detected in X1: mm = ", mm1, ", m = ", 
                m1, call. = FALSE)
    dummy <- simpleRDA2(Y, X2, SS.Y, mm2)
    bc.ua <- dummy$Rsquare
    m2 <- dummy$m
    if (m2 != mm2) 
        warning("collinearity detected in X2: mm = ", mm2, ", m = ", 
                m2, call. = FALSE)
    mm3 <- mm1 + mm2
    dummy <- simpleRDA2(Y, cbind(X1, X2), SS.Y, mm3)
    abc.ua <- dummy$Rsquare
    m3 <- dummy$m
    if (m3 != mm3) 
        warning("collinearity detected in cbind(X1,X2): mm = ", 
                mm3, ", m = ", m3, call. = FALSE)
    if ((m1 + m2) > m3) 
        bigwarning <- c("X1, X2")
    else bigwarning <- NULL
    ab <- RsquareAdj(ab.ua, n, m1)
    bc <- RsquareAdj(bc.ua, n, m2)
    abc <- RsquareAdj(abc.ua, n, m3)
    Df <- c(m1, m2, m3)
    fract <- data.frame(Df = Df,
                        R.squared = c(ab.ua, bc.ua, abc.ua),
                        Adj.R.squared = c(ab, bc, abc),
                        Testable = rep(TRUE, 3) & Df)
    rownames(fract) <- c("[a+b] = X1", "[b+c] = X2", "[a+b+c] = X1+X2")
    b <- ab + bc - abc
    Df <- c(m3-m2, 0, m3-m1, NA)
    indfract <- data.frame(Df = Df, R.squared = rep(NA, 4),
                           Adj.R.squared = c(ab - b, b, bc - b, 1 - abc),
                           Testable = c(TRUE, FALSE, TRUE, FALSE) & Df)
    rownames(indfract) <- c("[a] = X1|X2", "[b]", "[c] = X2|X1", 
                            "[d] = Residuals")
    out <- list(SS.Y = SS.Y, fract = fract, indfract = indfract, 
                nsets = 2, bigwarning = bigwarning, n = n1)
    class(out) <- "varpart234"
    out
}

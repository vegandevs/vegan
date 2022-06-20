`varpart2` <-
    function (Y, X1, X2, chisquare, permat)
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
    n <- nrow(Y)
    n1 <- nrow(X1)
    n2 <- nrow(X2)
    mm1 <- ncol(X1)
    mm2 <- ncol(X2)
    if (n1 != n)
        stop("Y and X1 do not have the same number of rows")
    if (n2 != n)
        stop("Y and X2 do not have the same number of rows")
    X1 <- scale(X1, center = TRUE, scale = FALSE)
    X2 <- scale(X2, center = TRUE, scale = FALSE)
    dummy <- simpleRDA2(Y, X1, SS.Y, permat)
    ab.ua <- dummy$Rsquare
    ab <- dummy$RsquareAdj
    m1 <- dummy$m
    if (m1 != mm1)
        collinwarn("X1", mm1, m1)
    dummy <- simpleRDA2(Y, X2, SS.Y, permat)
    bc.ua <- dummy$Rsquare
    bc <- dummy$RsquareAdj
    m2 <- dummy$m
    if (m2 != mm2)
        collinwarn("X2", mm2, m2)
    mm3 <- mm1 + mm2
    dummy <- simpleRDA2(Y, cbind(X1, X2), SS.Y, permat)
    abc.ua <- dummy$Rsquare
    abc <- dummy$RsquareAdj
    m3 <- dummy$m
    if (m3 != mm3)
        collinwarn("cbind(X1,X2)", mm3, m3)
    if ((m1 + m2) > m3)
        bigwarning <- c("X1, X2")
    else bigwarning <- NULL
    Df <- c(m1, m2, m3)
    ## labels and variable names are inconsistent below: the var names
    ## are from the old model where independent fractions were [a] and
    ## [c] and shared fraction was [b], but this was made consistent
    ## with other models where independent fractions are [a], [b],
    ## [c], [d] and listed before shared fractions.
    fract <- data.frame(Df = Df,
                        R.squared = c(ab.ua, bc.ua, abc.ua),
                        Adj.R.squared = c(ab, bc, abc),
                        Testable = rep(TRUE, 3) & Df)
    rownames(fract) <- c("[a+c] = X1", "[b+c] = X2", "[a+b+c] = X1+X2")
    b <- ab + bc - abc
    Df <- c(m3-m2, m3-m1, 0, NA)
    indfract <- data.frame(Df = Df, R.squared = rep(NA, 4),
                           Adj.R.squared = c(ab - b, bc - b, b, 1 - abc),
                           Testable = c(TRUE, TRUE, FALSE, FALSE) & Df)
    rownames(indfract) <- c("[a] = X1|X2", "[b] = X2|X1", "[c]",
                            "[d] = Residuals")
    out <- list(SS.Y = SS.Y, fract = fract, indfract = indfract,
                nsets = 2, bigwarning = bigwarning, n = n1)
    class(out) <- "varpart234"
    out
}

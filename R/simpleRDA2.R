"simpleRDA2" <-
function (Y, X, SS.Y, ...)
{
    Q <- qr(X, tol=1e-6)
    Yfit.X <- qr.fitted(Q, Y)
    SS <- sum(Yfit.X^2)
    if (missing(SS.Y)) SS.Y <- sum(Y^2)
    Rsquare <- SS/SS.Y
    list(Rsquare = Rsquare, m = Q$rank)
}


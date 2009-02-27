permat.tsdiag <-
function(x, type = "bray") {
    tsVec <- ts(summary(x)[[type]])
    ARmod <- arima(tsVec, order = c(1, 0, 0))
    ARmod$call <- match.call()
    ARresid <- ARmod$residuals
    BOX.tsVec <- Box.test(tsVec, lag = 1, type = "Ljung-Box")
    BOX.tsVec$data.name <- switch(type,
        "bray" = "Bray-Curtis dissimilarities",
        "chisq" = "Chi-squared values")
    BOX.ARresid <- Box.test(ARresid, lag = 1, type = "Ljung-Box")
    BOX.ARresid$data.name <- "AR residuals"
    x$perm <- NULL
    out <- list(call=match.call(), x=x, ts=tsVec,
        arima=ARmod, box.ts=BOX.tsVec, box.resid=BOX.ARresid)
    class(out) <- "permat.tsdiag"
    out
}

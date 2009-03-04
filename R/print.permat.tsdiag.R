print.diagnose.permat <-
function(x, ...) {
    cat("First order autoregressive model\n")
    print(x$arima)
    cat("\nTest of independence\n")
    print(x$box.ts)
    print(x$box.resid)
}

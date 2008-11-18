"summary.procrustes" <-
function (object, digits = getOption("digits"), ...) 
{
    ans <- object[c("call", "ss")]
    n <- nrow(object$Yrot)
    k <- ncol(object$Yrot)
    ans$resid <- residuals(object)
    rmse <- sqrt(object$ss/n)
    ans$n <- n
    ans$k <- k
    ans$rmse <- rmse
    ans$rotation <- object$rotation
    ans$translation <- object$translation
    ans$scale <- object$scale
    ans$digits <- digits
    class(ans) <- "summary.procrustes"
    ans
}


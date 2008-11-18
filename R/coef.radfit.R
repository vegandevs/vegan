"coef.radfit" <-
    function (object, ...) 
{
    out <- sapply(object$models, function(x) if (length(coef(x)) < 
                                                 3) 
                  c(coef(x), rep(NA, 3 - length(coef(x))))
    else coef(x))
    out <- t(out)
    colnames(out) <- paste("par", 1:3, sep = "")
    out
}

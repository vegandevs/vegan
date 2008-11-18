"summary.humpfit" <-
    function (object, ...) 
{
    dispersion <- if (any(object$family$family == c("binomial", 
                          "poisson"))) 
        1
    else sum(object$residuals^2)/object$df.residual
    p <- coef(object)
    se <- sqrt(dispersion * diag(solve(object$nlm$hessian)))
    est <- cbind(p, se)
    colnames(est) <- c("Estimate", "Std. Error")
    covmat <- solve(object$nlm$hessian)
    dg <- sqrt(diag(covmat))
    cormat <- covmat/outer(dg, dg)
    colnames(cormat) <- names(p)
    rownames(cormat) <- names(p)
    aic <- AIC(object)
    bic <- AIC(object, k = log(length(object$y)))
    out <- list(est = est, aic = aic, bic = bic, family = family(object)$family, 
                deviance = deviance(object), df.residual = df.residual(object), 
                dispersion = dispersion, correlation = cormat,
                cov.unscaled = covmat,  iter = object$nlm$iterations,
                code = object$nlm$code)
    class(out) <- "summary.humpfit"
    out
}

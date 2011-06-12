`anova.prc` <-
    function(object, ...)
{
    ## Refit prc() as an rda() for anova.cca()
    Z <- qr.X(object$pCCA$QR)
    X <- qr.X(object$CCA$QR)
    ## Get the name of the original response variable
    fla <- as.character(object$call$response)
    fla <- as.formula(paste(fla, "~ X + Condition(Z)"))
    mod <- rda(fla)
    anova(mod, ...)
}

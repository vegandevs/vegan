### 'working' will be Gower's G = -GowerDblcen(dis^2)/2
`fitted.dbrda` <-
    function (object, model = c("CCA", "CA", "pCCA"),
              type = c("response", "working"), ...)
{
    ZAP <- sqrt(.Machine$double.eps)
    type <- match.arg(type)
    model <- match.arg(model)
    if (is.null(object[[model]]))
        stop(gettextf("component '%s' does not exist", model))
    D <- ordiYbar(object, model)
    if (type == "response") {
        ## revert Gower double centring
        de <- diag(D)
        D <- -2 * D + outer(de, de, "+")
        ## we may have tiny negative zeros: zero them, but let large
        ## negative values be and give NaN in sqrt (with a warning)
        D[abs(D) < ZAP] <- 0
        if (!object$sqrt.dist)
            D <- sqrt(D)
        D <- D * object$adjust
        D <- as.dist(D)
        ## we do not remove Lingoes or Cailliez adjustment: this
        ## typically gives too many negative distances as unadjusted D
        ## often has zero-values
    }
    D
}

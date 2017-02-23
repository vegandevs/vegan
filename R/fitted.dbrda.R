### 'working' will be Gower's G = -GowerDblcen(dis^2)/2
`fitted.dbrda` <-
    function (object, model = c("CCA", "CA", "pCCA"),
              type = c("response", "working"), ...)
{
    ZAP <- sqrt(.Machine$double.eps)
    type <- match.arg(type)
    model <- match.arg(model)
    if (is.null(object[[model]]))
        stop("component ", model, " does not exist")
    D <- ordiYbar(object, model)
    if (type == "response") {
        ## revert Gower double centring
        de <- diag(D)
        D <- -2 * D + outer(de, de, "+")
        ## we may have tiny negative zeros
        if (any(D < 0)) {
            D[abs(D) < ZAP] <- 0
            if (any(D) < 0)
                warning("some squared dissimilarities are negative")
        }
        D <- sqrt(D)
        D <- D * object$adjust
        D <- as.dist(D)
    }
    D
}

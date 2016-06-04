### 'working' will be Gower's G = -GowerDblcen(dis^2)/2
`fitted.dbrda` <-
    function (object, model = c("CCA", "CA", "pCCA"),
              type = c("response", "working"), ...) 
{
    type <- match.arg(type)
    model <- match.arg(model)
    if (object$adjust == 1)
        const <- nobs(object) - 1
    else
        const <- 1
    if (is.null(object[[model]]))
        stop("component ", model, " does not exist")
    if (type == "working") {
        if (model == "pCCA")
            G <- object$pCCA$Fit
        else
            G <- object[[model]]$G
        if (model == "CCA") {
            H <- tcrossprod(
                qr.Q(object$CCA$QR)[, seq_len(object$CCA$QR$rank),
                                    drop=FALSE])
            G <- H %*% G %*% H
        }
        out <- G
    }
    if (type == "response") {
        if (model == "pCCA")
            stop("type = 'response' is unavailable for 'pCCA'")
        eig <- object[[model]]$eig
        U <- object[[model]]$u
        U <- sweep(U, 2, sqrt(eig[eig>0]), "*")
        D <- dist(U)
        ## remove additive constant
        if (!is.null(object$ac)) {
            if (object$add == "lingoes")
                D <- sqrt(D^2 - 2 * object$ac)
            else if (object$add == "cailliez")
                D <- D - object$ac
            else stop("unknown Euclidifying adjustment")
        }
        ## remove negative distances in imaginary space
        if (any(eig < 0)) {
            U <- object[[model]]$imaginary.u
            U <- sweep(U, 2, sqrt(abs(eig[eig<0])), "*")
            D <- sqrt(D^2 - dist(U)^2)
        }
        ## undo internal sqrt.dist
        if (object$sqrt.dist)
            D <- D^2
        out <- D * sqrt(const)
    }
    out
}

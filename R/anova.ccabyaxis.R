`anova.ccabyaxis` <-
    function (object, cutoff = 1,  ...) 
{
    cutoff <- cutoff + sqrt(.Machine$double.eps)
    rnk <- object$CCA$rank
    if (!max(rnk, 0)) 
        stop("Needs a constrained ordination")
    if (is.null(object$terms)) 
        stop("Analysis is only possible for models fitted using formula")
    ## Handle missing values in scores, both "omit" and "exclude" to
    ## match dims with data.
    if (!is.null(object$na.action)) {
        u <- stats:::napredict.exclude(object$na.action, object$CCA$u)
    } else {
        u <- object$CCA$u
    }
    ## Get conditions
    if (!is.null(object$pCCA)) {
        CondMat <- qr.X(object$pCCA$QR)
        ## deweight if CCA
        if (!inherits(object, "rda"))
            CondMat <- sweep(CondMat, 1, sqrt(object$rowsum), "/")
    }
    else
        CondMat <- NULL
    ## pad with NA rows if there is a subset
    if (!is.null(object$subset)) {
        lc <- matrix(NA, nrow=length(object$subset),
                     ncol = NCOL(u))
        lc[object$subset,]  <- u
        if (!is.null(CondMat)) {
            tmp <- matrix(NA, nrow=length(object$subset),
                          ncol = NCOL(CondMat))
            tmp[object$subset,] <- CondMat
            CondMat <- tmp
        }
        object$call$subset <- object$subset
    } else {
        lc <- u
    }
    lc <- as.data.frame(lc)
    axnam <- colnames(lc)
    df <- c(rep(1, rnk), object$CA$rank)
    chi <- c(object$CCA$eig, Residual = object$CA$tot.chi)
    Fval <- c(chi[1:rnk]/df[1:rnk]/chi[rnk+1]*df[rnk+1], NA)
    nperm <- c(numeric(rnk), NA)
    Pval <- rep(NA, rnk+1)
    out <- data.frame(df, chi, Fval, nperm, Pval)
    environment(object$terms) <- environment()
    fla <- paste(". ~ ", axnam[1], "+ Condition(",
                 paste(axnam[-1], collapse="+"),")")
    if (!is.null(CondMat)) {
        fla <- paste(fla, " + Condition(CondMat)")
        lc$CondMat <- CondMat
    }
    fla <- update(formula(object), fla)
    sol <- anova(update(object, fla, data=lc),  ...)
    out[c(1, rnk + 1), ] <- sol
    seed <- attr(sol, "Random.seed")
    attr(out, "names") <- attr(sol, "names")
    .call <- pasteCall(object$call, "Model:")
    attr(out, "heading") <- sub(" \n","", .call)
    attr(out, "Random.seed") <- seed
    bigseed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    bigperm <- out$N.Perm[1]
    if (rnk > 1) {
        for (.ITRM in 2:rnk) {
            fla <- paste(".~", axnam[.ITRM], "+Condition(",
                         paste(axnam[-(.ITRM)], collapse="+"),")")
            if (!is.null(CondMat))
                fla <- paste(fla, "+ Condition(CondMat)")
            fla <- update(formula(object),  fla) 
            sol <- update(object, fla, data = lc)
            assign(".Random.seed", seed, envir = .GlobalEnv)
            out[.ITRM, ] <- as.matrix(anova(sol, ...))[1,]
            if (out[.ITRM, "N.Perm"] > bigperm) {
                bigperm <- out[.ITRM, "N.Perm"]
                bigseed <- get(".Random.seed", envir = .GlobalEnv, 
                  inherits = FALSE)
            }
            if (out[.ITRM, "Pr(>F)"] > cutoff)
                break
        }
    }
    assign(".Random.seed", bigseed, envir = .GlobalEnv)
    class(out) <- c("anova.cca", "anova", "data.frame")
    out
}

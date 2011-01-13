"anova.ccabyaxis" <-
function (object, cutoff = 1,  ...) 
{
    cutoff <- cutoff + sqrt(.Machine$double.eps)
    rnk <- object$CCA$rank
    if (!max(rnk, 0)) 
        stop("Needs a constrained ordination")
    if (is.null(object$terms)) 
        stop("Analysis is only possible for models fitted using formula")
    lc <- object$CCA$u
    ## Handle missing values in scores, both "omit" and "exclude" to
    ## match dims with data.
    if (!is.null(object$na.action)) {
        lc <- stats:::napredict.exclude(object$na.action, lc)
    }
    axnam <- colnames(lc)
    df <- c(rep(1, rnk), object$CA$rank)
    chi <- c(object$CCA$eig, Residual = object$CA$tot.chi)
    Fval <- c(chi[1:rnk]/df[1:rnk]/chi[rnk+1]*df[rnk+1], NA)
    nperm <- c(numeric(rnk), NA)
    Pval <- rep(NA, rnk+1)
    out <- data.frame(df, chi, Fval, nperm, Pval)
    environment(object$terms) <- environment()
    fla <- update(formula(object), . ~ lc[,1] + Condition(lc[,-1]))
    sol <- anova(update(object, fla),  ...)
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
            fla <- update(formula(object),  .~ lc[, .ITRM] + Condition(lc[,-(.ITRM)]) )
            sol <- update(object, fla)
            assign(".Random.seed", seed, envir = .GlobalEnv)
            out[.ITRM, ] <- as.matrix(anova(sol, ...))[1, 
                ]
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

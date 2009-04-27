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
    newdata <- cbind(lc, eval(as.list(object$call)$data))
    axnam <- colnames(lc)
    df <- c(rep(1, rnk), object$CA$rank)
    chi <- c(object$CCA$eig, Residual = object$CA$tot.chi)
    Fval <- c(chi[1:rnk]/df[1:rnk]/chi[rnk+1]*df[rnk+1], NA)
    nperm <- c(numeric(rnk), NA)
    Pval <- rep(NA, rnk+1)
    out <- data.frame(df, chi, Fval, nperm, Pval)
    sol <- anova(object, first = TRUE, ...)
    out[c(1, rnk + 1), ] <- sol
    seed <- attr(sol, "Random.seed")
    attr(out, "names") <- attr(sol, "names")
    attr(out, "heading") <- attr(sol, "heading")
    attr(out, "Random.seed") <- seed
    bigseed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    bigperm <- out$N.Perm[1]
    environment(object$terms) <- environment()
    if (rnk > 1) {
        for (.ITRM in 2:rnk) {
            zz <- paste(paste("Condition(", axnam[1:(.ITRM - 1)], 
                ")"), collapse = "+")
            fla <- update(formula(object), paste(". ~ . +", zz))
            sol <- update(object, fla, data = newdata)
            assign(".Random.seed", seed, envir = .GlobalEnv)
            out[.ITRM, ] <- as.matrix(anova(sol, first = TRUE, ...))[1, 
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

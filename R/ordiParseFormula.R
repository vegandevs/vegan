## Pass thru Conditions or partial terms
Condition <- function(x) x

ordiParseFormula <- function(formula, data, xlev = NULL, na.action = na.fail,
                             subset = NULL, X)
{
    ## get terms
    trms0 <- terms(formula, specials = "Condition", data = data)
    ## Evaluate response data and delete from terms
    if (missing(X))
        Y <- as.matrix(eval(trms0[[2]], environment(formula),
                            enclos=.GlobalEnv))
    else
        Y <- as.matrix(X)
    trms <- delete.response(trms0)
    ## see if there are any Conditions partialled out
    pterm <- attr(trms, "specials")$Condition
    ## evaluate data
    mf <- model.frame(trms, data = data,
                      na.action = na.pass, # evaluate only after subset
                      drop.unused.levels = TRUE,
                      xlev = xlev)
    ## subset
    if (!is.null(subset)) {
        subset <- eval(subset,
                       if (is.data.frame(data)) cbind(data, Y)
                       else as.data.frame(Y),
                       parent.frame(2))
        ## subset will drop 'terms'
        att <- terms(mf)
        mf <- subset(mf, subset, drop = FALSE)
        attr(mf, "terms") <- att
        Y <- subset(Y, subset, drop = FALSE)
    }
    ## na.action
    if (any(is.na(mf))) {
        mf <- na.action(mf)
        nas <- attr(mf, "na.action")
        excluded <- Y[nas,, drop=FALSE]
        Y <- Y[-nas,, drop=FALSE]
    } else
        excluded <- NULL
    ## Separate X and Z
    trmlab <- attr(trms, "term.labels")
    ## empty model: unconstrained ordination
    if (trms[[2]] == "1" || trms[[2]] == "0") {
        X <- NULL
        Z <- NULL
    } else if (length(pterm) == 0 || is.null(pterm)) {
        X <- model.matrix(reformulate(trmlab), mf)
        Z <- NULL
    } else {
        X <- model.matrix(reformulate(trmlab[-pterm]), mf)
        Z <- model.matrix(reformulate(trmlab[pterm]), mf)[,-1,drop=FALSE]
    }
    ## Intercept is removed, but this also removes assign argument
    if (NROW(X) > 0) {
        assign <- attr(X, "assign")
        assign <- assign[assign > 0]
        X <- X[, -1, drop = FALSE]
        attr(X, "assign") <- assign
    }
    list(X = Y, Y = X, Z = Z, terms = trms0, terms.expand = trms0,
         modelframe = mf, subset = subset, na.action = attr(mf, "na.action"),
         excluded = excluded)
}

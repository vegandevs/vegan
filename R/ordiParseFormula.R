## Pass thru Conditions or partial terms
Condition <- function(x) x

ordiParseFormula <- function(formula, data, xlev = NULL, na.action = na.fail,
                             subset = NULL, X)
{
    ## Not yet implemented arguments
    if (!missing(xlev))
        .NotYetUsed("xlev")
    if (!missing(X))
        .NotYetUsed("X")
    ## get terms
    trms0 <- terms(formula, specials = "Condition", data = data)
    ## Evaluate response data and delete from terms
    Y <- as.matrix(eval(trms0[[2]]))
    trms <- delete.response(trms0)
    ## see if there are any Conditions partialled out
    pterm <- attr(trms, "specials")$Condition
    ## evaluate data
    mf <- model.frame(trms, data = data,
                      na.action = na.pass, # evaluate only after subset
                      drop.unused.levels = TRUE)
    ## subset
    if (!is.null(subset)) {
        subset <- eval(subset,
                       if (is.data.frame(data)) cbind(data, Y)
                       else as.data.frame(Y),
                       parent.frame(2))
        mf <- subset(mf, subset, drop = FALSE)
        Y <- subset(Y, subset, drop = FALSE)
    }
    ## na.action
    if (any(is.na(mf))) {
        mf <- na.action(mf)
        nas <- attr(mf, "na.action")
        Y <- Y[-nas,, drop=FALSE]
    }
    ## Separate X and Z
    trmlab <- attr(trms, "term.labels")
    if (length(pterm) == 0 || is.null(pterm)) {
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
         modelframe = mf)
}

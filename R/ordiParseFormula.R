ordiParseFormula <- function(formula, data, xlev = NULL, na.action = na.fail,
                             subset = NULL, X)
{
    ## Not yet implemented arguments
    if (!missing(xlev))
        .NotYetUsed("xlev")
    if (!missing(na.action))
        .NotYetUsed("na.action", error = FALSE)
    if (!missing(subset))
        .NotYetUsed("subset", error = FALSE)
    if (!missing(X))
        .NotYetUsed("X")
    ## Pass thru Conditions or partial terms
    Condition <- function(x) x
    ## get terms
    trms0 <- terms(formula, specials = "Condition", data = data)
    ## Evaluate response data and delete from terms
    Y <- as.matrix(eval(trms0[[2]]))
    trms <- delete.response(trms0)
    ## see if there are any Conditions partialled out
    pterm <- attr(trms, "specials")$Condition
    ## evaluate data
    mf <- model.frame(trms, data = data)
    ## Separate X and Y
    trmlab <- attr(trms, "term.labels")
    if (length(pterm) == 0 || is.null(pterm)) {
        X <- model.matrix(reformulate(trmlab), mf)[,-1,drop=FALSE]
        Z <- NULL
    } else {
        X <- model.matrix(reformulate(trmlab[-pterm]), mf)[,-1,drop=FALSE]
        Z <- model.matrix(reformulate(trmlab[pterm]), mf)[,-1,drop=FALSE]
    }
    list(X = Y, Y = X, Z = Z, terms = trms0, terms.expand = trms0,
         modelframe = mf)
}

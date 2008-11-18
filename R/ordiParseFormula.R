"ordiParseFormula" <-
function (formula, data, xlev = NULL, envdepth = 2) 
{
    Terms <- terms(formula, "Condition", data = data)
    flapart <- fla <- formula <- formula(Terms, width.cutoff = 500)
    specdata <- formula[[2]]
    X <- eval.parent(specdata, n = envdepth)
    X <- as.matrix(X)
    indPartial <- attr(Terms, "specials")$Condition
    mf <- Z <- NULL
    if (!is.null(indPartial)) {
        partterm <- attr(Terms, "variables")[1 + indPartial]
        Pterm <- sapply(partterm, function(x) deparse(x[[2]], width.cutoff=500))
        Pterm <- paste(Pterm, collapse = "+")
        P.formula <- as.formula(paste("~", Pterm), env = environment(formula))
        zlev <- xlev[names(xlev) %in% Pterm]
        mf <- model.frame(P.formula, data, na.action = na.fail, 
            xlev = zlev)
        Z <- model.matrix(P.formula, mf)
        if (any(colnames(Z) == "(Intercept)")) {
            xint <- which(colnames(Z) == "(Intercept)")
            Z <- Z[, -xint, drop = FALSE]
        }
        partterm <- sapply(partterm, function(x) deparse(x, width.cutoff=500))
        formula <- update(formula, paste(".~.-", paste(partterm, 
            collapse = "-")))
        flapart <- update(formula, paste(". ~ . +", Pterm))
    }
    formula[[2]] <- NULL
    if (formula[[2]] == "1" || formula[[2]] == "0") 
        Y <- NULL
    else {
        if (exists("Pterm")) 
            xlev <- xlev[!(names(xlev) %in% Pterm)]
        mf <- model.frame(formula, data, na.action = na.fail, 
            xlev = xlev)
        Y <- model.matrix(formula, mf)
        if (any(colnames(Y) == "(Intercept)")) {
            xint <- which(colnames(Y) == "(Intercept)")
            Y <- Y[, -xint, drop = FALSE]
        }
    }
    rownames(X) <- rownames(X, do.NULL = FALSE)
    colnames(X) <- colnames(X, do.NULL = FALSE)
    if (!is.null(Y)) {
        rownames(Y) <- rownames(Y, do.NULL = FALSE)
        colnames(Y) <- colnames(Y, do.NULL = FALSE)
    }
    if (!is.null(Z)) {
        rownames(Z) <- rownames(Z, do.NULL = FALSE)
        colnames(Z) <- colnames(Z, do.NULL = FALSE)
    }
    list(X = X, Y = Y, Z = Z, terms = terms(fla, width.cutoff = 500), 
        terms.expand = terms(flapart, width.cutoff = 500), modelframe = mf)
}

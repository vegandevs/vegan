`ordiParseFormula` <-
function (formula, data, xlev = NULL, na.action = na.fail,
          subset = NULL, X)
{
    if (missing(data))
        data <- environment(formula)
    fla <- formula
    ## distance-based methods (capscale, dbrda) evaluate specdata (LHS
    ## in formula) within their code, and the handling in
    ## ordiParseFormula is redundand and can be expensive
    if (missing(X)) {
        specdata <- formula[[2]]
        X <- eval(specdata, environment(formula), enclos=globalenv())
    }
    X <- as.matrix(X)
    formula[[2]] <- NULL
    Terms <- terms(formula, "Condition", data = data)
    indPartial <- attr(Terms, "specials")$Condition
    mf <- get_all_vars(formula, data)
    zmf <- ymf <- NULL
    Y <- Z <- NULL
    if (!is.null(indPartial)) {
        zmf <- mf[, indPartial, drop=FALSE]
        if (ncol(zmf) < ncol(mf))
            ymf <- mf[, -indPartial, drop=FALSE]
    } else {
        ymf <- mf
    }
    if (formula[[2]] == "1" || formula[[2]] == "0")
        Y <- NULL
    else {
        if (!is.null(zmf))
            xlev <- xlev[names(xlev) %in% names(ymf)]
    }
    ## Combine condition an constrain data frames
    if (!is.null(zmf)) {
        ncond <- NCOL(zmf)
        if (!is.null(ymf))
            mf <- cbind(zmf, ymf)
        else
            mf <- zmf
    } else {
        ncond <- 0
        mf <- ymf
    }
    ## Select a subset of data and species
    if (!is.null(subset)) {
        subset <- eval(subset,
                       if (inherits(data, "data.frame")) cbind(data, X)
                       else as.data.frame(X),
                       parent.frame(2))
        X <- X[subset, , drop = FALSE]
        if (NROW(mf) > 0)
            mf <- mf[subset, , drop = FALSE]
    }
    ## Get na.action attribute, remove NA and drop unused levels
    if (NCOL(mf) > 0) {
        mf <- model.frame(formula(mf), mf, xlev = xlev,
                          na.action = na.action, drop.unused.levels = TRUE)
        nas <- attr(mf, "na.action")
        ## Check if there are one-level factors after subset and na.action
        for (i in 1:ncol(mf))
            if (is.factor(mf[[i]]) && length(levels(mf[[i]])) <= 1)
                levels(mf[[i]]) <- c(levels(mf[[i]]), ".ThisVarHasOnly1Level")
    } else {
        nas <- NULL
    }
    ## Check and remove NA in dependent data
    if (!is.null(nas)) {
        excluded <- X[nas, , drop = FALSE]
        X <- X[-nas,, drop=FALSE]
    } else {
        excluded <-  NULL
    }
    if (ncond > 0) {
        Z <- model.matrix(reformulate(names(zmf)), zmf)
        if (any(colnames(Z) == "(Intercept)"))
            Z <- Z[, -which(colnames(Z) == "(Intercept)"), drop = FALSE]
    }
    if (NCOL(ymf) > 0) {
        Y <- model.matrix(reformulate(names(ymf)), ymf)
        ## save assign attribute
        assign <- attr(Y, "assign")
        assign <- assign[assign > 0]
        if (any(colnames(Y) == "(Intercept)"))
            Y <- Y[, -which(colnames(Y) == "(Intercept)"), drop = FALSE]
        if (NCOL(Y) == 0)
            Y <- NULL
        else
            attr(Y, "assign") <- assign
    }
    X <- as.matrix(X)
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
    if (NCOL(mf) > 0)
        Trms <- terms(reformulate(names(mf)), width.cutoff = 500)
    else {
        Trms <- terms(fla)
        mf <- NULL
    }
    list(X = X, Y = Y, Z = Z, terms = fla,
         terms.expand = Trms, modelframe = mf,
         subset = subset, na.action = nas, excluded = excluded)
}

"ordiParseFormula" <-
function (formula, data, xlev = NULL, envdepth = 2, na.action = na.fail,
          subset = NULL) 
{
    if (missing(data))
        data <- environment(formula)
    Terms <- terms(formula, "Condition", data = data)
    flapart <- fla <- formula <- formula(Terms, width.cutoff = 500)
    specdata <- formula[[2]]
    X <- eval.parent(specdata, n = envdepth)
    indPartial <- attr(Terms, "specials")$Condition
    zmf <- ymf <- Y <- Z <- NULL
    formula[[2]] <- NULL
    if (!is.null(indPartial)) {
        partterm <- attr(Terms, "variables")[1 + indPartial]
        Pterm <- sapply(partterm, function(x) deparse(x[[2]], width.cutoff=500))
        Pterm <- paste(Pterm, collapse = "+")
        P.formula <- as.formula(paste("~", Pterm), env = environment(formula))
        zlev <- xlev[names(xlev) %in% Pterm]
        zmf <- model.frame(P.formula, data, na.action = na.pass, 
            xlev = zlev)
        partterm <- sapply(partterm, function(x) deparse(x, width.cutoff=500))
        formula <- update(formula, paste("~.-", paste(partterm, 
            collapse = "-")))
        flapart <- update(formula, paste(" ~ . +", Pterm))
    }
    if (formula[[2]] == "1" || formula[[2]] == "0") 
        Y <- NULL
    else {
        if (exists("Pterm")) 
            xlev <- xlev[!(names(xlev) %in% Pterm)]
        ymf <- model.frame(formula, data, na.action = na.pass, 
            xlev = xlev)
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
                       if (inherits(data, "data.frame")) cbind(data, X) else X,
                       parent.frame())
        X <- X[subset, , drop = FALSE]
        if (NROW(mf) > 0)
            mf <- mf[subset, , drop = FALSE]
    }
    ## Get na.action attribute, remove NA and drop unused levels
    if (NROW(mf) > 0) {
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
        Z <- model.matrix(P.formula, mf)
        if (any(colnames(Z) == "(Intercept)"))
            Z <- Z[, -which(colnames(Z) == "(Intercept)"), drop = FALSE]
    }
    if (NROW(mf) > 0) {
        Y <- model.matrix(formula, mf)
        if (any(colnames(Y) == "(Intercept)"))
            Y <- Y[, -which(colnames(Y) == "(Intercept)"), drop = FALSE]
        if (NCOL(Y) == 0)
            Y <- NULL
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
    list(X = X, Y = Y, Z = Z, terms = terms(fla, width.cutoff = 500), 
         terms.expand = terms(flapart, width.cutoff = 500), modelframe = mf,
         na.action = nas, excluded = excluded)
}

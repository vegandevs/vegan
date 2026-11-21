`ordiParseFormula` <-
function (formula, data, xlev = NULL, na.action = na.fail,
          subset = NULL, X)
{
    if (missing(data))
        data <- environment(formula)
    Terms <- terms(formula, "Condition", data = data)
    flapart <- fla <- formula <- formula(Terms, width.cutoff = 500)
    ## distance-based methods (capscale, dbrda) evaluate specdata (LHS
    ## in formula) within their code, and the handling in
    ## ordiParseFormula is redundant and can be expensive
    if (missing(X)) {
        specdata <- formula[[2]]
        X <- eval(specdata, environment(formula), enclos=globalenv())
    }
    X <- as.matrix(X, rownames.force = TRUE)
    indPartial <- attr(Terms, "specials")$Condition
    zmf <- ymf <- Y <- Z <- NULL
    formula[[2]] <- NULL
    if (!is.null(indPartial)) {
        partterm <- attr(Terms, "variables")[1 + indPartial]
        ## paste() needed to combine >500 character deparsed Conditions
        Pterm <- sapply(partterm, function(x)
            paste(deparse(x[[2]], width.cutoff=500, backtick = TRUE), collapse=" "))
        Pterm <- paste(Pterm, collapse = "+")
        P.formula <- as.formula(paste("~", Pterm), env = environment(formula))
        zlev <- xlev[names(xlev) %in% Pterm]
        zmf <- if (inherits(data, "environment"))
            eval(substitute(
                model.frame(P.formula, na.action = na.pass, xlev = zlev)),
                 envir = data, enclos = .GlobalEnv)
        else
            model.frame(P.formula, data, na.action = na.pass, xlev = zlev)
        partterm <- sapply(partterm, function(x)
            paste(deparse(x, width.cutoff=500, backtick = TRUE), collapse = " "))
        formula <- update(formula, paste("~.-", paste(partterm,
            collapse = "-")))
        flapart <- update(formula, paste(" ~ . +", Pterm))
    }
    if (formula[[2]] == "1" || formula[[2]] == "0")
        Y <- NULL
    else {
        if (exists("Pterm"))
            xlev <- xlev[!(names(xlev) %in% Pterm)]

        ymf <- if (inherits(data, "environment"))
            eval(substitute(
                model.frame(formula, na.action = na.pass, xlev = xlev)),
                 envir=data, enclos=.GlobalEnv)
        else
            model.frame(formula, data, na.action = na.pass, xlev = xlev)
    }
    ## Record prediction calls before subsetting so their transformation
    ## settings come from the same data as the fitted model matrix.
    constraint.terms <- if (!is.null(ymf)) ordiPredictTerms(ymf, data) else NULL
    ## Combine the condition and constraint model frames
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
    if (NROW(mf) > 0) {
        ## We need to re-construct model.frame and its formula, and
        ## for this names with functions must be back-quoted
        ## (`poly(A1, 2)`, e.g.) with the current scoping to find the
        ## variable (A1) later. This happened with formula(mf) in R <
        ## 3.6.0, but now we must be explicit.
        mf <- model.frame(as.formula(paste("~",
              paste(sapply(names(mf), "as.name"), collapse="+"))),
              mf, xlev = xlev,
              na.action = na.action, drop.unused.levels = TRUE)
        nas <- attr(mf, "na.action")
        ## Check if there are one-level factors after subset and na.action
        for (i in 1:ncol(mf))
            if (is.factor(mf[[i]]) && nlevels(mf[[i]]) <= 1)
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
        ## save assign attribute
        assign <- attr(Y, "assign")
        contrasts <- attr(Y, "contrasts")
        assign <- assign[assign > 0]
        if (any(colnames(Y) == "(Intercept)"))
            Y <- Y[, -which(colnames(Y) == "(Intercept)"), drop = FALSE]
        if (NCOL(Y) == 0)
            Y <- NULL
        else {
            attr(Y, "assign") <- assign
            attr(Y, "contrasts") <- contrasts
        }
    }
    X <- as.matrix(X) # if X was dist
    rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "") # if X was matrix
    colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "spe")
    list(X = X, Y = Y, Z = Z, terms = terms(fla, width.cutoff = 500),
         terms.expand = terms(flapart, width.cutoff = 500), modelframe = mf,
         constraint.terms = constraint.terms,
         subset = subset, na.action = nas, excluded = excluded)
}

## model.frame records top-level makepredictcall results in predvars.
## Extend these to fitted transformations inside standard wrappers.
ordiPredictTerms <- function(mf, data)
{
    Terms <- attr(mf, "terms")
    vars <- attr(Terms, "predvars")
    if (!is.null(vars)) {
        original <- attr(Terms, "variables")
        for (i in seq_along(mf)) {
            j <- i + 1L
            vars[j] <- list(ordiPredictCall(vars[[j]], data, environment(Terms),
                                            value = mf[[i]],
                                            original = original[[j]]))
        }
        attr(Terms, "predvars") <- vars
    }
    Terms
}

ordiPredictCall <- function(expr, data, env, value = NULL, original = expr)
{
    if (missing(expr))
        return(quote(expr = ))
    if (!is.call(expr))
        return(expr)
    head <- expr[[1L]]
    name <- if (is.symbol(head)) as.character(head)
            else if (is.call(head) &&
                     is.symbol(head[[1L]]) &&
                     as.character(head[[1L]]) %in% c("::", ":::"))
                as.character(head[[3L]]) else ""
    wrappers <- c("I", "[", "(", "+", "-", "*", "/", "^")
    if (!name %in% c(wrappers, "scale", "scale.default", "poly", "ns", "bs"))
        return(expr)
    ## Look up the function as R did while fitting, skipping data columns
    ## with the same name. Leave arbitrary calls to makepredictcall methods:
    ## their arguments may be quoted or evaluated in another environment.
    fun <- if (is.symbol(head))
        eval(substitute(base::get(NAME, mode = "function"), list(NAME = name)),
             envir = data, enclos = env)
    else eval(head, envir = data, enclos = env)
    is.scale <- identical(fun, base::scale) || identical(fun, base::scale.default)
    is.poly <- identical(fun, stats::poly)
    is.spline <- !is.null(fun) && name %in% c("ns", "bs") &&
                 identical(environmentName(environment(fun)), "splines")
    if (!is.scale && !is.poly && !is.spline &&
        !(name %in% wrappers && identical(fun, get(name, baseenv()))))
        return(expr)
    if (is.scale) {
        ## makepredictcall() can append named settings to positional arguments.
        ## Start with the original call so matching does not use them twice.
        expr <- match.call(fun, if (is.null(value)) expr else original)
    }
    for (i in seq_along(expr)[-1L]) {
        ## I() wraps its result without changing the transformation attributes.
        inner <- if (identical(fun, base::I)) value else NULL
        expr[i] <- list(ordiPredictCall(expr[[i]], data, env, value = inner))
    }
    if (!is.scale && !is.poly && !is.spline)
        return(expr)
    ## Keep fitted calls already supplied by model.frame(). Only their
    ## nested arguments may still need fitted settings of their own.
    if (is.poly) {
        if (is.null(value) && !"coefs" %in% names(expr))
            value <- eval(expr, envir = data, enclos = env)
        if (!is.null(coefs <- attr(value, "coefs")))
            expr$coefs <- coefs
        else if ("coefs" %in% names(expr))
            coefs <- eval(expr$coefs, envir = data, enclos = env)
        ## Univariate coefficients can also come from a one-column matrix,
        ## e.g. poly(scale(x), 2). Use the vector prediction path without
        ## evaluating the input again to inspect its dimensions.
        if (is.list(coefs) && all(c("alpha", "norm2") %in% names(coefs))) {
            expr <- match.call(stats::poly, expr)
            expr$x <- as.call(list(quote(base::as.numeric), expr$x))
        }
        return(expr)
    }
    if (is.spline && all(c("knots", "Boundary.knots") %in% names(expr)))
        return(expr)
    if (is.null(value))
        value <- eval(expr, envir = data, enclos = env)
    if (is.scale) {
        if (!is.null(center <- attr(value, "scaled:center")))
            expr$center <- center
        if (!is.null(scale <- attr(value, "scaled:scale")))
            expr$scale <- scale
        return(expr)
    }
    stats::makepredictcall(value, expr)
}

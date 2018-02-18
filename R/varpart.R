`varpart` <-
    function (Y, X, ..., data, chisquare = FALSE, transfo, scale = FALSE,
              add = FALSE, sqrt.dist = FALSE, permutations = 1000)
{
    if (missing(data))
        data <- parent.frame()
    X <- list(X, ...)
    if ((length(X) < 2 || length(X) > 4))
        stop("needs two to four explanatory tables")
    ## transfo and scale can be used only with non-distance data
    if (inherits(Y, "dist")) {
        inert <- attr(Y, "method")
        inert <- paste(paste0(toupper(substring(inert, 1, 1)),
                              substring(inert, 2)), "distance")
        ## sqrt of distances?
        if (sqrt.dist)
            Y <- sqrt(Y)
        else
            inert <- paste("squared", inert)
        ## additive constant to euclidify distances?
        if (is.logical(add) && isTRUE(add))
            add <- "lingoes"
        if (is.character(add)) {
            add <- match.arg(add, c("lingoes", "cailliez"))
            if (add == "lingoes") {
                ac <- addLingoes(as.matrix(Y))
                Y <- sqrt(Y^2 + 2 * ac)
            } else if (add == "cailliez") {
                ac <- addCailliez(as.matrix(Y))
                Y <- Y + ac
            }
            if (ac > sqrt(.Machine$double.eps))
                inert <- paste(paste0(toupper(substring(add, 1, 1)),
                                    substring(add, 2)), "adjusted", inert)
        }
        RDA <- "dbRDA"
        if(!missing(transfo) || !missing(scale))
            message("arguments 'transfo' and 'scale' are ignored with distances")
    } else if (chisquare) {
        inert = "Chi-square"
        RDA = "CCA"
        permutations = getPermuteMatrix(permutations, nrow(Y))
    } else {
        inert <- "variance"
        RDA <- "RDA"
        if (!missing(transfo)) {
            Y <- decostand(Y, transfo)
            transfo <- attr(Y, "decostand")
        }
        if (!missing(transfo) && (is.null(dim(Y)) || ncol(Y) == 1))
            warning("transformations are probably meaningless with a single variable")
        if (scale && !missing(transfo))
            warning("Y should not be both transformed and scaled (standardized)")
        Y <- scale(Y, center = TRUE, scale = scale)
    }
    Sets <- list()
    for (i in seq_along(X)) {
        if (is.data.frame(X[[i]]) || is.factor(X[[i]])) {
            ## factor variable or a data.frame (possibly with factors)
            mf <- as.data.frame(X[[i]])
            mf <- model.matrix(~ ., mf)
            Sets[[i]] <- mf[,-1, drop = FALSE] # remove intercept
        } else if (inherits(X[[i]], "formula")) {
            ## Formula interface
            mf <- model.frame(X[[i]], data, na.action = na.fail,
                              drop.unused.levels = TRUE)
            trms <- attr(mf, "terms")
            Sets[[i]] <- model.matrix(trms, mf)
            if (any(colnames(Sets[[i]]) == "(Intercept)")) {
                xint <- which(colnames(Sets[[i]]) == "(Intercept)")
                Sets[[i]] <- (Sets[[i]])[, -xint, drop = FALSE]
            }
        }
        else Sets[[i]] <- as.matrix(X[[i]])
        Sets[[i]] <- scale(Sets[[i]], center = TRUE, scale = TRUE)
    }
    out <- list()
    out$part <- switch(length(Sets), NULL,
                       varpart2(Y, Sets[[1]], Sets[[2]],
                                chisquare, permutations),
                       varpart3(Y, Sets[[1]], Sets[[2]], Sets[[3]],
                                chisquare, permutations),
                       varpart4(Y, Sets[[1]], Sets[[2]], Sets[[3]], Sets[[4]],
                                chisquare, permutations))
    if (inherits(Y, "dist")) {
        out$part$ordination <- "dbrda"
    } else {
        if (chisquare)
            out$part$ordination <- "cca"
        else {
            out$part$ordination <- "rda"
        }
    }
    if(RDA == "RDA") {
        out$scale <- scale
        if (!missing(transfo))
            out$transfo <- transfo
    } else {
        if (scale || !missing(transfo))
            message("arguments 'scale' and 'transfo' ignored: valid only in RDA")
    }
    out$inert <- inert
    out$RDA <- RDA
    out$call <- match.call()
    mx <- rep(" ", length(X))
    for (i in seq_along(X)) mx[i] <- deparse(out$call[[i+2]], width.cutoff = 500)
    out$tables <- mx
    class(out) <- c("varpart", class(out))
    out
}

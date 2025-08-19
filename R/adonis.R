`adonis2` <-
    function(formula, data, permutations = 999, method = "bray",
             sqrt.dist = FALSE, add = FALSE, by = NULL,
             parallel = getOption("mc.cores"), na.action = na.fail,
             strata = NULL, ...)
{
    ## handle missing data
    if (missing(data))
        data <- parent.frame()
    else
        data <- eval(match.call()$data, environment(formula),
                     enclos = .GlobalEnv)
    formula <- formula(terms(formula, data = data))
    ## we accept only by = "terms", "margin", "onedf" or NULL
    if (!is.null(by))
        by <- match.arg(by, c("terms", "margin", "onedf"))
    ## evaluate lhs
    lhs <- eval(formula[[2]], envir = environment(formula),
                enclos = globalenv())
    ## Take care that input lhs are dissimilarities
    if ((is.matrix(lhs) || is.data.frame(lhs)) &&
        isSymmetric(unname(as.matrix(lhs))))
        lhs <- as.dist(lhs)
    if (!inherits(lhs, "dist"))
        lhs <- vegdist(as.matrix(lhs), method=method, ...)
    ## evaluate terms like in dbrda
    d <- ordiParseFormula(formula = formula,
                          data = data,
                          na.action = na.action,
                          subset = NULL,
                          X = lhs) # do not re-evaluate lhs
    if (is.null(d$Y))
        stop("needs explanatory variables on the right-hand-side")
    lhs <- d$X
    ## handle missing data
    if (!is.null(d$na.action))
        lhs <- lhs[-d$na.action, -d$na.action, drop = FALSE]
    ## adjust distances if requested
    if (sqrt.dist)
        lhs <- sqrt(lhs)
    if (is.logical(add) && add)
        add <- "lingoes"
    if (is.character(add)) {
        add <- match.arg(add, c("lingoes", "cailliez"))
        if (add == "lingoes") {
            ac <- addLingoes(as.matrix(lhs))
            lhs <- sqrt(lhs^2 + 2 * ac)
        }
        else if (add == "cailliez") {
            ac <- addCailliez(as.matrix(lhs))
            lhs <- lhs + ac
        }
    }
    sol <- adonis0(lhs, d$Y, d$Z)
    sol$formula <- match.call()
    sol$terms <- d$terms
    sol$terminfo <- ordiTerminfo(d, data)
    ## handle permutations
    perm <- getPermuteMatrix(permutations, NROW(lhs), strata = strata)
    out <- anova(sol, permutations = perm, by = by,
                 parallel = parallel)
    ## attributes will be lost when adding a new column
    att <- attributes(out)
    ## add traditional adonis output on R2
    out <- rbind(out, "Total" = c(nobs(sol)-1, sol$tot.chi, NA, NA))
    out <- cbind(out[,1:2], "R2" = out[,2]/sol$tot.chi, out[,3:4])
    ## Fix output header to show the adonis2() call instead of adonis0()
    att$heading[2] <- deparse(match.call(), width.cutoff = 500L)
    att$names <- names(out)
    att$row.names <- rownames(out)
    attributes(out) <- att
    out
}

## adonis0 does the same as ordConstrained, except ordination
`adonis0` <-
    function(lhs, X, Z)
{
    ## G is -dmat/2 centred
    G <- initDBRDA(lhs)
    ## output object
    sol <- list(Ybar = G,
                tot.chi = sum(diag(G)),
                adjust = 1,
                method = "adonis")
    ## partial with Condition
    if (!is.null(Z) && ncol(Z)) {
        Z <- scale(Z, scale = FALSE)
        QZ <- qr(Z)
        G <- qr.resid(QZ, G)
        G <- qr.resid(QZ, t(G))
        pCCA <- list(rank = QZ$rank,
                     tot.chi = rowSums(qr.fitted(QZ, G)),
                     QR = QZ)
    } else {
        pCCA <- NULL
    }
    ## always have constraints X if we got here
    if (!is.null(Z))
        X <- cbind(Z, X)
    X <- scale(X, scale = FALSE)
    qrhs <- qr(X)
    Gfit <- qr.fitted(qrhs, G)
    Gfit <- qr.fitted(qrhs, t(Gfit))
    Gres <- qr.resid(qrhs, G)
    Gres <- qr.resid(qrhs, t(Gres))
    ## collect data for the fit
    rank <- qrhs$rank
    if (!is.null(Z))
        rank <- rank - QZ$rank
    if(rank > 0)
        CCA <- list(rank = rank,
                    qrank = rank,
                    tot.chi = sum(diag(Gfit)),
                    QR = qrhs)
    else
        CCA <- NULL # empty model
    ## collect data for the residuals
    CA <- list(rank = nrow(lhs) - max(qrhs$rank, 0) - 1,
               u = matrix(0, nrow = nrow(lhs)),
               tot.chi = sum(diag(Gres)))
    ## all together
    if (!is.null(pCCA))
        sol$pCCA <- pCCA
    sol$CCA <- CCA
    sol$CA <- CA
    class(sol) <- c("adonis2", "dbrda", "rda", "cca")
    sol
}

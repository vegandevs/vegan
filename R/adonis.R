`adonis` <-
    function(formula, data=NULL, permutations=999, method="bray", strata=NULL,
             contr.unordered="contr.sum", contr.ordered="contr.poly",
             parallel = getOption("mc.cores"), ...)
{
    ## formula is model formula such as Y ~ A + B*C where Y is a data
    ## frame or a matrix, and A, B, and C may be factors or continuous
    ## variables.  data is the data frame from which A, B, and C would
    ## be drawn.
    TOL <- 1e-7
    Terms <- terms(formula, data = data)
    lhs <- formula[[2]]
    lhs <- eval(lhs, data, parent.frame()) # to force evaluation
    formula[[2]] <- NULL                # to remove the lhs
    rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE) # to get the data frame of rhs
    op.c <- options()$contrasts
    options( contrasts=c(contr.unordered, contr.ordered) )
    rhs <- model.matrix(formula, rhs.frame) # and finally the model.matrix
    options(contrasts=op.c)
    grps <- attr(rhs, "assign")
    qrhs <- qr(rhs)
    ## Take care of aliased variables and pivoting in rhs
    rhs <- rhs[, qrhs$pivot, drop=FALSE]
    rhs <- rhs[, 1:qrhs$rank, drop=FALSE]
    grps <- grps[qrhs$pivot][1:qrhs$rank]
    u.grps <- unique(grps)
    nterms <- length(u.grps) - 1
    if (nterms < 1)
        stop("right-hand-side of formula has no usable terms")
    H.s <- lapply(2:length(u.grps),
                  function(j) {Xj <- rhs[, grps %in% u.grps[1:j] ]
                               qrX <- qr(Xj, tol=TOL)
                               Q <- qr.Q(qrX)
                               tcrossprod(Q[,1:qrX$rank])
                           })
    if (inherits(lhs, "dist")) {
        if (any(lhs < -TOL))
            stop("dissimilarities must be non-negative")
        dmat <- as.matrix(lhs^2)
    }
    else {
        dist.lhs <- as.matrix(vegdist(lhs, method=method, ...))
        dmat <- dist.lhs^2
    }
    n <- nrow(dmat)
    ## G is -dmat/2 centred by rows
    G <- -sweep(dmat, 1, rowMeans(dmat))/2
    SS.Exp.comb <- sapply(H.s, function(hat) sum( G * t(hat)) )
    SS.Exp.each <- c(SS.Exp.comb - c(0,SS.Exp.comb[-nterms]) )
    H.snterm <- H.s[[nterms]]
    ## t(I - H.snterm) is needed several times and we calculate it
    ## here
    tIH.snterm <- t(diag(n)-H.snterm)
    if (length(H.s) > 1)
        for (i in length(H.s):2)
            H.s[[i]] <- H.s[[i]] - H.s[[i-1]]
    SS.Res <- sum( G * tIH.snterm)
    df.Exp <- sapply(u.grps[-1], function(i) sum(grps==i) )
    df.Res <- n - qrhs$rank
    ## Get coefficients both for the species (if possible) and sites
    if (inherits(lhs, "dist")) {
        beta.sites <- qr.coef(qrhs, as.matrix(lhs))
        beta.spp <-  NULL
    } else {
        beta.sites <- qr.coef(qrhs, dist.lhs)
        beta.spp <-  qr.coef(qrhs, as.matrix(lhs))
    }
    colnames(beta.spp) <- colnames(lhs)
    colnames(beta.sites) <- rownames(lhs)
    F.Mod <- (SS.Exp.each/df.Exp) / (SS.Res/df.Res)

    f.test <- function(tH, G, df.Exp, df.Res, tIH.snterm) {
      ## HERE I TRY CHANGING t(H)  TO tH, and
      ## t(I - H.snterm) to tIH.snterm, so that we don't have
      ## to do those calculations for EACH iteration.
      ## This is the function we have to do for EACH permutation.
      ## G is an n x n centered distance matrix
      ## H is the hat matrix from the design (X)
      ## note that for R, * is element-wise multiplication,
      ## whereas %*% is matrix multiplication.
        (sum(G * tH)/df.Exp) /
          (sum(G * tIH.snterm)/df.Res) }

 ### Old f.test
    ### f.test <- function(H, G, I, df.Exp, df.Res, H.snterm){
    ##    (sum( G * t(H) )/df.Exp) /
      ##    (sum( G * t(I-H.snterm) )/df.Res) }

    SS.perms <- function(H, G, I){
        c(SS.Exp.p = sum( G * t(H) ),
          S.Res.p=sum( G * t(I-H) )
          ) }

    ## Permutations
    p <- getPermuteMatrix(permutations, n, strata = strata)
    permutations <- nrow(p)
    if (permutations) {
        tH.s <- lapply(H.s, t)
        ## Apply permutations for each term
        ## This is the new f.test (2011-06-15) that uses fewer arguments
        ## Set first parallel processing for all terms
        if (is.null(parallel))
            parallel <- 1
        hasClus <- inherits(parallel, "cluster")
        isParal <- hasClus || parallel > 1
        isMulticore <- .Platform$OS.type == "unix" && !hasClus
        if (isParal && !isMulticore && !hasClus) {
            parallel <- makeCluster(parallel)
        }
        if (isParal) {
            if (isMulticore) {
                f.perms <-
                    sapply(1:nterms, function(i)
                           unlist(mclapply(1:permutations, function(j)
                                           f.test(tH.s[[i]], G[p[j,], p[j,]],
                                                  df.Exp[i], df.Res, tIH.snterm),
                                           mc.cores = parallel)))
            } else {
                f.perms <-
                    sapply(1:nterms, function(i)
                           parSapply(parallel, 1:permutations, function(j)
                                     f.test(tH.s[[i]], G[p[j,], p[j,]],
                                            df.Exp[i], df.Res, tIH.snterm)))
            }
        } else {
            f.perms <-
                sapply(1:nterms, function(i) 
                       sapply(1:permutations, function(j) 
                              f.test(tH.s[[i]], G[p[j,], p[j,]],
                                     df.Exp[i], df.Res, tIH.snterm)))
        }
        ## Close socket cluster if created here
        if (isParal && !isMulticore && !hasClus)
            stopCluster(parallel)
        ## Round to avoid arbitrary P-values with tied data
        f.perms <- round(f.perms, 12)
        P <- (rowSums(t(f.perms) >= F.Mod)+1)/(permutations+1)
    } else { # no permutations
        f.perms <- P <- rep(NA, nterms)
    }
    F.Mod <- round(F.Mod, 12)
    SumsOfSqs = c(SS.Exp.each, SS.Res, sum(SS.Exp.each) + SS.Res)
    tab <- data.frame(Df = c(df.Exp, df.Res, n-1),
                      SumsOfSqs = SumsOfSqs,
                      MeanSqs = c(SS.Exp.each/df.Exp, SS.Res/df.Res, NA),
                      F.Model = c(F.Mod, NA,NA),
                      R2 = SumsOfSqs/SumsOfSqs[length(SumsOfSqs)],
                      P = c(P, NA, NA))
    rownames(tab) <- c(attr(attr(rhs.frame, "terms"), "term.labels")[u.grps],
                       "Residuals", "Total")
    colnames(tab)[ncol(tab)] <- "Pr(>F)"
    attr(tab, "heading") <- c(howHead(attr(p, "control")),
        "Terms added sequentially (first to last)\n")
    class(tab) <- c("anova", class(tab))
    out <- list(aov.tab = tab, call = match.call(),
                coefficients = beta.spp, coef.sites = beta.sites,
                f.perms = f.perms, model.matrix = rhs, terms = Terms)
    class(out) <- "adonis"
    out
}

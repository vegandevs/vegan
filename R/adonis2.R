`adonis2` <-
    function(formula, data=NULL, method="bray", ...)
{
    ## First we collect info for the uppermost level of the analysed
    ## object
    Trms <- terms(delete.response(formula), data = data)
    sol <- list(call = match.call(),
                method = "adonis",
                formula = formula(Trms),
                terms = Trms,
                terminfo = list(terms = Trms))
    sol$call$formula <- formula(Trms)
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
    ##op.c <- options()$contrasts
    ##options( contrasts=c(contr.unordered, contr.ordered) )
    rhs <- model.matrix(formula, rhs.frame) # and finally the model.matrix
    ##grps <- attr(rhs, "assign")
    rhs <- rhs[,-1, drop=FALSE] # remove the (Intercept) to get rank right
    ##grps <- grps[-1]
    rhs <- scale(rhs, scale = FALSE, center = TRUE) # center
    ##options(contrasts=op.c)
    qrhs <- qr(rhs)
    ## handle dissimilarities
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
    ## G is -dmat/2 centred
    G <- -GowerDblcen(dmat)/2
    ## preliminaries are over: start working
    Gfit <- qr.fitted(qrhs, G)
    Gres <- qr.resid(qrhs, G)
    ## collect data for the fit
    if(!is.null(qrhs$rank) && qrhs$rank > 0) 
        CCA <- list(rank = qrhs$rank,
                    qrank = qrhs$rank,
                    tot.chi = sum(diag(Gfit)),
                    QR = qrhs,
                    G = G)
    else
        CCA <- NULL # empty model
    ## collect data for the residuals
    CA <- list(rank = n - max(qrhs$rank, 0) - 1,
               u = matrix(0, nrow=n),
               tot.chi = sum(diag(Gres)),
               Xbar = Gres)
    ## all together
    sol$tot.chi <- sum(diag(G))
    sol$CCA <- CCA
    sol$CA <- CA
    class(sol) <- c("adonis2", "capscale", "rda", "cca")
    sol
}

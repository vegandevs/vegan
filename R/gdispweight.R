### Clarke's dispweight is based on the hypothesis that count data
### should follow Poisson distribution, and species overdispersed to
### the Poisson should be downweighted. The basic model assesses the
### expected values of species and their overdispersion wrt to class
### means for a single factor and then estimates the significance of
### the overdispersion using individual-based simulation within these
### same classes. Function gdispweight generalizes this by allowing a
### formula that specifies any fitted model, but estimates the
### significance of the overdispersion analytically from Pearson
### residuals.

`gdispweight` <-
    function(formula, data, plimit = 0.05)
{
    ## We do not handle missing values (yet?)
    op <- options(na.action = "na.fail")
    on.exit(op)
    ## extract response data
    comm <- eval(formula[[2]])
    ## extract rhs
    if (missing(data))
        data <- environment(formula)
    x <- model.matrix(delete.response(terms(formula, data = data)),
                      data = data)
    ## Quasi-Poisson
    family <- quasipoisson()
    V <- family$variance
    ## fit models to all species separately and extract results
    comm <- as.data.frame(comm)
    mods <- lapply(comm, function(y) glm.fit(x, y, family = family))
    y <- sapply(mods, '[[', "y")
    mu <- sapply(mods, fitted)
    wts <- sapply(mods, '[[',  "prior.weights")
    res <- (y-mu) * sqrt(wts) / sqrt(V(mu))
    df <- sapply(mods, df.residual)
    ## the same stats as in Clarke's original, but parametrically
    stat <- colSums(res^2)
    p <- pchisq(stat, df, lower.tail = FALSE)
    dhat <- stat/df
    w <- ifelse(p < plimit, 1/dhat, 1)
    ## do not upweight underdispersed species
    w <- ifelse(w > 1, 1, w)
    ## done
    comm <- sweep(comm, 2, w, "*")
    class(comm) <- c("dispweight", class(comm))
    attr(comm, "D") <- dhat
    attr(comm, "df") <- df
    attr(comm, "p") <- p
    attr(comm, "weights") <- w
    attr(comm, "nsimul") <- NA
    attr(comm, "nullmodel") <- NA
    comm
}

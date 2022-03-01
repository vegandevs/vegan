`centroids.cca` <-
    function(x, mf, wt)
{
    if (is.null(mf) || is.null(x))
        return(NULL)
    facts <- sapply(mf, is.factor) | sapply(mf, is.character)
    if (!any(facts))
        return(NULL)
    mf <- mf[, facts, drop = FALSE]
    ## Explicitly exclude NA as a level
    mf <- droplevels(mf, exclude = NA)
    if (missing(wt))
        wt <- rep(1, nrow(mf))
    ind <- seq_len(nrow(mf))
    workhorse <- function(x, wt)
        colSums(x * wt) / sum(wt)
    ## As NA not a level, centroids only for non-NA levels of each factor
    tmp <- lapply(mf, function(fct)
                  tapply(ind, fct, function(i) workhorse(x[i,, drop=FALSE], wt[i])))
    tmp <- lapply(tmp, function(z) sapply(z, rbind))
    pnam <- labels(tmp)
    out <- NULL
    if (ncol(x) == 1) {
        nm <- unlist(sapply(pnam,
                            function(nm) paste(nm, names(tmp[[nm]]), sep="")),
                     use.names=FALSE)
        out <- matrix(unlist(tmp), nrow=1, dimnames = list(NULL, nm))
    } else {
        for (i in seq_along(tmp)) {
            colnames(tmp[[i]]) <- paste(pnam[i], colnames(tmp[[i]]),
                                        sep = "")
            out <- cbind(out, tmp[[i]])
        }
    }
    out <- t(out)
    colnames(out) <- colnames(x)
    out
}

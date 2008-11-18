"centroids.cca" <-
    function (x, mf, wt) 
{
    mf <- mf[, unlist(lapply(mf, is.factor)), drop = FALSE]
    if (ncol(mf) == 0) 
        return(NULL)
    if (missing(wt)) 
        wt <- rep(1, nrow(mf))
    x <- sweep(x, 1, wt, "*")
    workhorse <- function(mf, x, wt) {
        sw <- tapply(wt, mf, sum)
        swx <- apply(x, 2, tapply, mf, sum)
        sweep(swx, 1, sw, "/")
    }
    tmp <- lapply(mf, workhorse, x, wt)
    pnam <- labels(tmp)
    out <- NULL
    for (i in 1:length(tmp)) {
        rownames(tmp[[i]]) <- paste(pnam[i], rownames(tmp[[i]]), 
                                    sep = "")
        out <- rbind(out, tmp[[i]])
    }
    out
}

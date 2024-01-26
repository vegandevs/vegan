### SIMPER: contributions of species on overall dissimilarity with
### emphasis on among-group dissimilarities. Generate full
### dissimilarity matrix first and then subsample.

`simper` <-
    function(comm, group, permutations = 999, parallel = 1, ...)
{
    ## parallel processing not yet implemented
    if (!missing(parallel))
        .NotYetUsed("parallel", error = FALSE)
    EPS <- sqrt(.Machine$double.eps)
    comm <- as.matrix(comm)
    ## take lower triangle without as.dist overhead
    tri <- outer(seq_len(nrow(comm)), seq_len(nrow(comm)), ">")
    ## Species contributions of differences needed for every species,
    ## but denominator is constant. Bray-Curtis is actually
    ## manhattan/(mean(rowsums)) and this is the way we collect data
    rs <- rowSums(comm)
    rs <- outer(rs, rs, "+")[tri]
    spcontr <- sapply(seq_len(ncol(comm)),
        function(i) as.vector(vegdist(comm[, i, drop = FALSE], "man")))
    ## Bray-Curtis
    spcontr <- sweep(spcontr, 1, rs, "/")
    colnames(spcontr) <- colnames(comm)
    outlist <- NULL
    ## Averages of species contributions
    ## Case 1: overall differences without grouping
    if (missing(group) || length(unique(group)) == 1) {
        nperm <- 0
        permat <- NULL
        average <- colMeans(spcontr)
        overall <- sum(average)
        sdi <- apply(spcontr, 2, sd)
        ord <- order(average, decreasing = TRUE)
        cusum <- cumsum(average[ord])/overall
        outlist[["total"]] <- list(species = colnames(comm),
                                   average = average, overall = overall,
                                   sd = sdi, ratio = average/sdi,
                                   ava = NULL, avb = NULL, ord = ord,
                                   cusum = cusum, p = NULL)
    } else {
        ## Case 2: two or more groups
        comp <- t(combn(as.character(unique(group)), 2))
        ## data averages by group (do we need these?)
        spavg <- apply(comm, 2, function(x) tapply(x, group, mean))
        ## function to match constrasts
        contrmatch <- function(X, Y, patt)
            X != Y & X %in% patt & Y %in% patt
        for (i in seq_len(nrow(comp))) {
            tmat <- outer(group, group, FUN=contrmatch, patt=comp[i,])
            take <- tmat[tri]
            average <- colMeans(spcontr[take,,drop=FALSE])
            overall <- sum(average)
            sdi <- apply(spcontr[take,,drop=FALSE], 2, sd)
            ratio <- average/sdi
            ord <- order(average, decreasing = TRUE)
            cusum <- cumsum(average[ord])/overall
            ava <- spavg[comp[i,1],]
            avb <- spavg[comp[i,2],]
            ## Permutation tests for average
            permat <- getPermuteMatrix(permutations, nrow(comm))
            nperm <- nrow(permat)
            if (nperm) {
                Pval <- rep(1, ncol(comm))
                for (k in seq_len(nperm)) {
                    take <- tmat[permat[k,],permat[k,]][tri]
                    Pval <- Pval + ((colMeans(spcontr[take,,drop = FALSE]) - EPS) >= average)
                }
                Pval <- Pval/(nperm+1)
                if (anyNA(ratio))
                    Pval[is.na(ratio)] <- NA
            } else {
                Pval <- NULL
            }
            ## output
            outlist[[paste(comp[i,], collapse="_")]] <-
                list(species = colnames(comm), average = average,
                     overall = overall, sd = sdi, ratio = ratio, ava = ava,
                     avb = avb, ord = ord, cusum = cusum, p = Pval)
        }
    }
    class(outlist) <- "simper"
    attr(outlist, "permutations") <- nperm
    attr(outlist, "control") <- attr(permat, "control")
    outlist
}
`print.simper` <-
    function(x, ...)
{
    cat("cumulative contributions of most influential species:\n\n")
    cusum <- lapply(x, function(z) z$cusum)
    spec <- lapply(x, function(z) z$species[z$ord])
    for (i in seq_along(cusum)) {
        names(cusum[[i]]) <- spec[[i]]
    }
    ## this probably fails with empty or identical groups that have 0/0 = NaN
    out <- lapply(cusum, function(z) z[seq_len(min(which(z >= 0.7)))])
    print(out)
    invisible(x)
}

`summary.simper` <-
    function(object, ordered = TRUE, digits = max(3, getOption("digits") - 3), ...)
{
    if (ordered) {
        out <- lapply(object, function(z)
            data.frame(cbind(average = z$average, sd = z$sd, ratio = z$ratio,
                       ava = z$ava, avb = z$avb)[z$ord, ]))
        cusum <- lapply(object, function(z) z$cusum)
        for(i in seq_along(out)) {
            out[[i]]$cumsum <- cusum[[i]]
            if(!is.null(object[[i]]$p)) {
                out[[i]]$p <- object[[i]]$p[object[[i]]$ord]
            }
        }
    }
    else {
        out <- lapply(object, function(z)
            data.frame(cbind(average = z$average, sd = z$sd, ratio = z$ratio,
                             ava = z$ava, avb = z$avb, p = z$p)))
    }
    attr(out, "digits") <- digits
    attr(out, "permutations") <- attr(object, "permutations")
    attr(out, "control") <- attr(object, "control")
    class(out) <- "summary.simper"
    out
}

`print.summary.simper`<-
    function(x, digits = attr(x, "digits"), ...)
{
    for (nm in names(x)) {
        cat("\nContrast:", nm, "\n\n")
        printCoefmat(x[[nm]], digits = digits, has.Pvalue = TRUE,
                     zap.ind = seq_len(ncol(x[[nm]])), ...)
    }
    if (!is.null(attr(x, "control")))
        cat(howHead(attr(x, "control")))
    invisible(x)
}


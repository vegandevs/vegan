`simper` <-
    function(comm, group, permutations = 0, trace = FALSE,
             parallel = getOption("mc.cores"), ...)
{
    EPS <- sqrt(.Machine$double.eps)
    if (any(rowSums(comm, na.rm = TRUE) == 0))
        warning("you have empty rows: results may be meaningless")
    pfun <- function(x, comm, comp, i, contrp) {
        groupp <- group[perm[x,]]
        ga <- comm[groupp == comp[i, 1], , drop = FALSE]
        gb <- comm[groupp == comp[i, 2], , drop = FALSE]
        n.a <- nrow(ga)
        n.b <- nrow(gb)
        for(j in seq_len(n.b)) {
            for(k in seq_len(n.a)) {
                mdp <- abs(ga[k, , drop = FALSE] - gb[j, , drop = FALSE])
                mep <- ga[k, , drop = FALSE] + gb[j, , drop = FALSE]
                contrp[(j-1)*n.a+k, ] <- mdp / sum(mep)
            }
        }
        colMeans(contrp)
    }
    comm <- as.matrix(comm)
    comp <- t(combn(unique(as.character(group)), 2))
    outlist <- NULL
    ## data parameters
    P <- ncol(comm)
    nobs <- nrow(comm)
    ## Make permutation matrix
    perm <- getPermuteMatrix(permutations, nobs, ...)
    ## check dims (especially if permutations was a matrix)
    if (ncol(perm) != nobs)
        stop(gettextf("'permutations' have %d columns, but data have %d rows",
                          ncol(perm), nobs))
    ## OK: take number of permutations
    nperm <- nrow(perm)
    if (nperm > 0)
        perm.contr <- matrix(nrow=P, ncol=nperm)
    ## Parallel processing ?
    if (is.null(parallel))
        parallel <- 1
    hasClus <- inherits(parallel, "cluster")
    isParal <- hasClus || parallel > 1
    isMulticore <- .Platform$OS.type == "unix" && !hasClus
    if (isParal && !isMulticore && !hasClus) {
        parallel <- makeCluster(parallel)
    }
    for (i in seq_len(nrow(comp))) {
        group.a <- comm[group == comp[i, 1], , drop = FALSE]
        group.b <- comm[group == comp[i, 2], , drop = FALSE]
        n.a <- nrow(group.a)
        n.b <- nrow(group.b)
        contr <- matrix(ncol = P, nrow = n.a * n.b)
        for (j in seq_len(n.b)) {
            for (k in seq_len(n.a)) {
                md <- abs(group.a[k, , drop = FALSE] - group.b[j, , drop = FALSE])
                me <- group.a[k, , drop = FALSE] + group.b[j, , drop = FALSE]
                contr[(j-1)*n.a+k, ] <- md / sum(me)
            }
        }
        average <- colMeans(contr)

        ## Apply permutations
        if(nperm > 0){
            if (trace)
                cat("Permuting", paste(comp[i,1], comp[i,2], sep = "_"), "\n")
            contrp <- matrix(ncol = P, nrow = n.a * n.b)

            if (isParal) {
                if (isMulticore){
                    perm.contr <- mclapply(seq_len(nperm), function(d)
                        pfun(d, comm, comp, i, contrp), mc.cores = parallel)
                    perm.contr <- do.call(cbind, perm.contr)
                } else {
                    perm.contr <- parSapply(parallel, seq_len(nperm), function(d)
                        pfun(d, comm, comp, i, contrp))
                }
            } else {
                perm.contr <- sapply(1:nperm, function(d)
                    pfun(d, comm, comp, i, contrp))
            }
            p <- (rowSums(apply(perm.contr, 2, function(x) x >= average - EPS)) + 1) / (nperm + 1)
        }
        else {
          p <- NULL
        }

        overall <- sum(average)
        sdi <- apply(contr, 2, sd)
        ratio <- average / sdi
        ava <- colMeans(group.a)
        avb <- colMeans(group.b)
        ord <- order(average, decreasing = TRUE)
        cusum <- cumsum(average[ord] / overall)
        out <- list(species = colnames(comm), average = average,
                    overall = overall, sd = sdi, ratio = ratio, ava = ava,
                    avb = avb, ord = ord, cusum = cusum, p = p)
        outlist[[paste(comp[i,1], "_", comp[i,2], sep = "")]] <- out
    }
    ## Close socket cluster if created here
    if (isParal && !isMulticore && !hasClus)
        stopCluster(parallel)
    attr(outlist, "permutations") <- nperm
    attr(outlist, "control") <- attr(perm, "control")
    class(outlist) <- "simper"
    outlist
}

### An alternative implementation. I expect this be much faster, in
### particular in permutation tests. Generate full dissimilarity
### matrix first and then subsample.

`simper2` <-
    function(comm, group, permutations = 0, ...)
{
    EPS <- sqrt(.Machine$double.eps)
    comm <- as.matrix(comm)
    ## Species contributions of differences needed for every species,
    ## but denominator is constant. Bray-Curtis is actually
    ## manhattan/(mean(rowsums)) and this is the way we collect data
    rs <- rowSums(comm)
    rs <- as.dist(outer(rs, rs, "+"))
    spcontr <- sapply(seq_len(ncol(comm)),
                      function(i) vegdist(comm[,i,drop=FALSE], "man"))
    ## Bray-Curtis
    spcontr <- sweep(spcontr, 1, rs, "/")
    colnames(spcontr) <- colnames(comm)
    outlist <- NULL
    ## Averages of species contributions
    ## Case 1: overall differences without grouping
    if (missing(group) || length(unique(group)) == 1) {
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
        ## take lower triangle without as.dist overhead
        tri <- outer(seq_along(group), seq_along(group), ">")
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
                    Pval <- Pval + ((colMeans(spcontr[take,]) - EPS) >= average)
                }
                Pval <- Pval/(nperm+1)
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


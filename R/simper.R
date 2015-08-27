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
            data.frame(contr = z$average, sd = z$sd, ratio = z$ratio,
                       av.a = z$ava, av.b = z$avb)[z$ord, ])
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
            data.frame(cbind(contr = z$average, sd = z$sd, 'contr/sd' = z$ratio,
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
    signif.stars <- getOption("show.signif.stars") && attr(x, "permutations") > 0
    starprint <- function(z) {
        if (signif.stars && any(z$p < 0.1)) {
            stars <- symnum(z$p, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                            symbols = c("***", "**", "*", ".", " "))
            z <- cbind(z, " " = format(stars))
        }
        z
    }
    out <- lapply(x, starprint)
    for (nm in names(out)) {
        cat("\nContrast:", nm, "\n\n")
        print(out[[nm]], digits = digits, ...)
    }
    if (signif.stars && any(sapply(x, function(z) z$p) < 0.1)) {
        leg <- attr(symnum(1, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                            symbols = c("***", "**", "*", ".", " ")), "legend")
        cat("---\nSignif. codes: ", leg, "\n")
    }
    if (!is.null(attr(x, "control")))
        cat(howHead(attr(x, "control")))
    invisible(x)
}


`bioenv.default` <-
function (comm, env, method = "spearman", index = "bray", upto = ncol(env), 
              trace = FALSE, partial = NULL, parallel = getOption("mc.cores"),
          ...) 
{
    if (is.null(partial)) {
        corfun <- function(dx, dy, dz, method) {
            cor(dx, dy, method=method)
        }
    } else {
        corfun <- function(dx, dy, dz, method) {
            rxy <- cor(dx, dy, method=method)
            rxz <- cor(dx, dz, method=method)
            ryz <- cor(dy, dz, method=method)
            (rxy - rxz*ryz)/sqrt(1-rxz*rxz)/sqrt(1-ryz*ryz)
        }
    }
    if (!is.null(partial))
        partpart <- deparse(substitute(partial))
    else
        partpart <- NULL
    if (!is.null(partial) && !inherits(partial, "dist"))
        partial <- dist(partial)
    if (!is.null(partial) && !pmatch(method, c("pearson", "spearman"), nomatch=FALSE))
        stop("method ", method, " invalid in partial bioenv")
    n <- ncol(env)
    ntake <- 2^n - 1
    ndone <- 0
    if (n > 8 || trace) {
        if (upto < n) 
            cat("Studying", nall <- sum(choose(n, 1:upto)), "of ")
        cat(ntake, "possible subsets (this may take time...)\n")
        flush.console()
    }
    x <- scale(env)
    best <- list()
    if (inherits(comm, "dist")) {
        comdis <- comm
        index <- attr(comdis, "method")
        if (is.null(index))
            index <- "unspecified"
    } else if (is.matrix(comm) && nrow(comm) == ncol(comm) &&
             isTRUE(all.equal(comm, t(comm)))) {
        comdis <- as.dist(comm)
        index <- "supplied square matrix"
    } else {
        comdis <- vegdist(comm, method = index)
    }
    ## Prepare for parallel processing
    if (is.null(parallel) && getRversion() >= "2.15.0")
        parallel <- get("default", envir = parallel:::.reg)
    if (is.null(parallel) || getRversion() < "2.14.0")
        parallel <- 1
    hasClus <- inherits(parallel, "cluster")
    isParal <- (hasClus || parallel > 1) && require(parallel)
    isMulticore <- .Platform$OS.type == "unix" && !hasClus
    if (isParal && !isMulticore && !hasClus) {
        parallel <- makeCluster(parallel)
    }
    ## get the number of clusters
    if (inherits(parallel, "cluster"))
        nclus <- length(parallel)
    else
        nclus <- parallel
    CLUSLIM <- 8
    ## The proper loop
    for (i in 1:upto) {
        if (trace) {
            nvar <- choose(n, i)
            cat("No. of variables ", i, ", No. of sets ", nvar, 
                "...", sep = "")
            flush.console()
        }
        sets <- t(combn(1:n, i))
        if (!is.matrix(sets)) 
            sets <- as.matrix(t(sets))
        if (isParal && nrow(sets) >= CLUSLIM*nclus) {
            if (isMulticore) {
                est <- unlist(mclapply(1:nrow(sets), function(j)
                                       corfun(comdis, dist(x[,sets[j, ]]), partial,
                                              method = method),
                                       mc.cores = parallel))
            } else {
                est <- parSapply(parallel, 1:nrow(sets), function(j)
                                  corfun(comdis, dist(x[,sets[j, ]]), partial,
                                         method = method))
            }
        } else {
            est <- sapply(1:nrow(sets), function(j) 
                          corfun(comdis, dist(x[,sets[j, ]]), partial,
                                 method = method))
        }
        best[[i]] <- list(best = sets[which.max(est), ], est = max(est))
        if (trace) {
            ndone <- ndone + nvar
            cat(" done (", round(100 * ndone/ntake, 1), "%)\n", 
                sep = "")
            flush.console()
        }
    }
    out <- list(names = colnames(env), method = method, index = index, 
                upto = upto, models = best, partial = partpart)
    out$call <- match.call()
    out$call[[1]] <- as.name("bioenv")
    class(out) <- "bioenv"
    out
}

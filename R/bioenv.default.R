`bioenv.default` <-
function (comm, env, method = "spearman", index = "bray", upto = ncol(env), 
          trace = FALSE, partial = NULL,
          metric = c("euclidean", "mahalanobis", "manhattan", "gower"),
          parallel = getOption("mc.cores"),
          ...) 
{
    metric <- match.arg(metric)
    method <- match.arg(method, eval(formals(cor)$method))
    if (any(sapply(env, is.factor)) && metric != "gower")
        stop("you have factors in 'env': only 'metric = \"gower\"' is allowed")
    if (is.null(partial)) {
        corfun <- function(dx, dy, dz, method, ...) {
            cor(dx, dy, method=method, ...)
        }
    } else {
        corfun <- function(dx, dy, dz, method, ...) {
            rxy <- cor(dx, dy, method=method, ...)
            rxz <- cor(dx, dz, method=method, ...)
            ryz <- cor(dy, dz, method=method, ...)
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
    upto <- min(upto, n)
    if (n > 8 || trace) {
        if (upto < n) 
            cat("Studying", nall <- sum(choose(n, 1:upto)), "of ")
        cat(ntake, "possible subsets (this may take time...)\n")
        flush.console()
    }
    ## Check metric and adapt data and distance function
    if (metric == "euclidean") {
        x <- scale(env, scale = TRUE)
        distfun <- function(x) dist(x)
    } else if (metric == "mahalanobis") {
        x <- as.matrix(scale(env, scale = FALSE))
        distfun <- function(x) dist(veganMahatrans(x))
    } else if (metric == "gower") {
        x <- env
        distfun <- function(x) daisy(x, metric = "gower")
    } else if (metric == "manhattan") {
        x <- decostand(env, "range")
        distfun <- function(x) dist(x, "manhattan")
    } else {
        stop("unknown metric")
    }
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
    if (is.null(parallel))
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
                                       corfun(comdis,
                                              distfun(x[,sets[j,],drop = FALSE]),
                                              partial, method = method, ...),
                                       mc.cores = parallel))
            } else {
                est <- parSapply(parallel, 1:nrow(sets), function(j)
                                  corfun(comdis, distfun(x[,sets[j,],drop = FALSE]),
                                         partial, method = method, ...))
            }
        } else {
            est <- sapply(1:nrow(sets), function(j) 
                          corfun(comdis, distfun(x[,sets[j,], drop=FALSE ]),
                                 partial, method = method, ...))
        }
        best[[i]] <- list(best = sets[which.max(est), ], est = max(est))
        if (trace) {
            ndone <- ndone + nvar
            cat(" done (", round(100 * ndone/ntake, 1), "%)\n", 
                sep = "")
            flush.console()
        }
    }
    whichbest <- which.max(lapply(best, function(tmp) tmp$est))
    out <- list(names = colnames(env), method = method, index = index,
                metric = metric, upto = upto, models = best,
                whichbest = whichbest,
                partial = partpart, x = x, distfun = distfun)
    out$call <- match.call()
    out$call[[1]] <- as.name("bioenv")
    class(out) <- "bioenv"
    out
}

## Function to extract the environmental distances used within
## bioenv. The default is to take the best model, but any model can be
## specified by its number.

`bioenvdist`  <-
    function(x, which = "best")
{
    ## any non-numeric argument is regarded as "best"
    if(!is.numeric(which))
        which <- x$whichbest
    x$distfun(x$x[, x$models[[which]]$best, drop = FALSE])
}



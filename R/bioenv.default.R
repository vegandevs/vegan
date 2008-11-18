`bioenv.default` <-
function (comm, env, method = "spearman", index = "bray", upto = ncol(env), 
              trace = FALSE, partial = NULL, ...) 
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
    comdis <- vegdist(comm, method = index)
    for (i in 1:upto) {
        if (trace) {
            nvar <- choose(n, i)
            cat("No. of variables ", i, ", No. of sets ", nvar, 
                "...", sep = "")
            flush.console()
        }
        sets <- ripley.subs(i, 1:n)
        if (!is.matrix(sets)) 
            sets <- as.matrix(t(sets))
        est <- numeric(nrow(sets))
        for (j in 1:nrow(sets)) est[j] <- corfun(comdis, dist(x[, 
                                                                sets[j, ]]), partial, method = method)
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

#' Permutation test for the area of convex hull or ellipse in ordination
#'
#' Finds if the area covered by a convex hull or fitted ellipse is
#' smaller than expected under null hypothesis using permutation test.
#'
#' @param ord 2-d ordination
#' @param factor defining groups
#' @param are of convex hull of or an ellipse
#' @param permutations: number, permutation matrix or a
#' \code{\link[permute]{how}} definition.
#' @param parallel parallel processing
#' @param \dots other parameters passed to area functions
#'
#' @author Jari Oksanen
`ordiareatest` <-
    function(ord, groups, area = c("hull", "ellipse"), permutations = 999,
             parallel = getOption("mc.cores"), ...)
{
    ## Function to find area
    area <- match.arg(area)
    areafun <- if (area == "hull") ordihull else ordiellipse
    areafun <- match.fun(areafun)
    ## Observed statistics
    obs <- summary(areafun(ord, groups, draw = "none", ...))["Area",]
    ## permutations
    pfun <- function(take, ...)
        summary(areafun(ord, groups[take], draw = "none", ...))["Area",]
    perm <- getPermuteMatrix(permutations, length(groups))
    nperm <- nrow(perm)
    if (is.null(parallel))
        parallel <- 1
    hasClus <- inherits(parallel, "cluster")
    if (hasClus || parallel > 1) {
        if(.Platform$OS.type == "unix" && !hasClus) {
            areas <- do.call(cbind,
                             mclapply(1:permutations,
                                      function(i, ...) pfun(perm[i,],...),
                                        mc.cores = parallel))
            } else {
                if (!hasClus) {
                    parallel <- makeCluster(parallel)
                }
                areas <- parApply(parallel, perm, MARGIN=1, pfun)
                if (!hasClus)
                    stopCluster(parallel)
            }
    } else {
        areas <- sapply(1:permutations, function(i, ...) pfun(perm[i,], ...))
    }
    signif <- (rowSums(areas <= obs) + 1)/(nperm + 1)
    out <- list("areas" = obs, "pvalues" = signif, "permutations" = areas,
                nperm = nperm, control = attr(perm, "control"), "kind" = area)
    class(out) <- "ordiareatest"
    out
}

### print method

`print.ordiareatest` <-
    function(x, ...)
{
    qu <- apply(x$permutations, 1, quantile, probs=c(0.05, 0.5))
    m <- cbind("Area" = x$areas, t(qu), "Pr(<sim)" = x$pvalues)
    cat("\n")
    cat(gettextf("Permutation test for the size of ordination %ss\nAlternative hypothesis: observed area is smaller than random %s\n\n", x$kind, x$kind))
    cat(howHead(x$control), "\n")
    printCoefmat(m, tst.ind=1:3)
    invisible(x)
}

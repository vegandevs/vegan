scoverage <- function(x, coverage.limit = 10) {
    x <- as.matrix(x)
    n <- rowSums(x)
    if (any(n <= 0))
        stop("negative or zero row totals not allowed")
    ## check if comm contains integer, especially for singletons
    if (any(x[x>0] < 1))
        warning("<1 non integer values detected: analysis might not be meaningful")
    if (abs(sum(x) - sum(as.integer(x))) > 10^-6)
        warning("non integer values detected")
    C1 <- 1 - rowSums(x==1) / n
    a <- decostand(x, "total")
    ifelse(x < coverage.limit, a * C1, a)
}

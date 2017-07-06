`diversity` <-
    function (x, index = "shannon", MARGIN = 1, base = exp(1))
{
    x <- drop(as.matrix(x))
    if (any(x < 0, na.rm = TRUE))
        stop("input data must be non-negative")
    INDICES <- c("shannon", "simpson", "invsimpson")
    index <- match.arg(index, INDICES)
    if (length(dim(x)) > 1) {
        total <- apply(x, MARGIN, sum)
        x <- sweep(x, MARGIN, total, "/")
    } else {
        x <- x/(total <- sum(x))
    }
    if (index == "shannon")
        x <- -x * log(x, base)
    else
        x <- x*x
    if (length(dim(x)) > 1)
        H <- apply(x, MARGIN, sum, na.rm = TRUE)
    else
        H <- sum(x, na.rm = TRUE)
    if (index == "simpson")
        H <- 1 - H
    else if (index == "invsimpson")
        H <- 1/H
    ## check NA in data
    if (any(NAS <- is.na(total)))
        H[NAS] <- NA
    H
}

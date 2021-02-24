`diversity` <-
    function (x, index = "shannon", groups, equalize.groups = FALSE,
              MARGIN = 1, base = exp(1))
{
    x <- drop(as.matrix(x))
    if (!is.numeric(x))
        stop("input data must be numeric")
    if (any(x < 0, na.rm = TRUE))
        stop("input data must be non-negative")
    ## sum communities for groups
    if (!missing(groups)) {
        if (MARGIN == 2)
            x <- t(x)
        if (length(groups) == 1) # total for all SU
            groups <- rep(groups, NROW(x))
        if (equalize.groups)
            x <- decostand(x, "total")
        x <- aggregate(x, list(groups), sum) # pool SUs by groups
        rownames(x) <- x[,1]
        x <- x[,-1, drop=FALSE]
        if (MARGIN == 2)
            x <- t(x)
    }
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

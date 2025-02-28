## evaluate user-defined dissimilarity function.

`designdist` <-
    function (x, method = "(A+B-2*J)/(A+B)",
              terms = c("binary", "quadratic", "minimum"),
              abcd = FALSE, alphagamma = FALSE, name, maxdist)
{
    terms <- match.arg(terms)
    if ((abcd || alphagamma) && terms != "binary")
        warning("perhaps terms should be 'binary' with 'abcd' or 'alphagamma'?")
    x <- as.matrix(x)
    ## only do numeric data for which "pa", minimum and quadratic make sense
    if (!(is.numeric(x) || is.logical(x)))
        stop("input data must be numeric")
    N <- nrow(x)
    P <- ncol(x)
    if (terms == "binary")
        x <- ifelse(x > 0, 1, 0)
    if (terms == "binary" || terms == "quadratic")
        XX <- tcrossprod(x)
    if (terms == "minimum")
        XX <- .Call(do_minterms, as.matrix(x))
    d <- diag(XX)
    A <- as.dist(outer(rep(1, N), d))
    B <- as.dist(outer(d, rep(1, N)))
    J <- as.dist(XX)
    ## 2x2 contingency table notation
    if (abcd) {
        a <- J
        b <- A - J
        c <- B - J
        d <- P - A - B + J
    }
    ## beta diversity notation
    if (alphagamma) {
        alpha <- (A + B)/2
        gamma <- A + B - J
        delta <- abs(A - B)/2
    }
    dis <- eval(parse(text = method))
    attributes(dis) <- attributes(J)
    attr(dis, "call") <- match.call()
    if (missing(name))
        attr(dis, "method") <- paste(terms, method)
    else attr(dis, "method") <- name
    if (!missing(maxdist)) {
        if (!is.na(maxdist) && any(dis > maxdist)) {
            warning("'maxdist' was lower than some distances: setting to NA")
            maxdist <- NA
        }
        attr(dis, "maxdist") <- maxdist
    }
    dis
}

### Variant of designdist that calculates dissimilarities between two
### data sets. Returns a rectangular matrix where rows are the rows of
### matrix 'x', and columns are the rows of matrix 'y'.

`designdist2` <-
    function (x, y, method = "(A+B-2*J)/(A+B)",
              terms = c("binary", "quadratic", "minimum"),
              abcd = FALSE, alphagamma = FALSE, name, maxdist)
{
    terms <- match.arg(terms)
    if ((abcd || alphagamma) && terms != "binary")
        warning("perhaps terms should be 'binary' with 'abcd' or 'alphagamma'?")
    x <- as.matrix(x)
    y <- as.matrix(y)
    ## only do numeric data for which "pa", minimum and quadratic make sense
    if (!(is.numeric(x) || is.logical(x)))
        stop("input data must be numeric")
    Nx <- nrow(x)
    Ny <- nrow(y)
    P <- ncol(x)
    if (ncol(x) != P)
        stop("number of columns do not match in 'x' and 'y'")
    if (terms == "binary") {
        x <- ifelse(x > 0, 1, 0)
        y <- ifelse(y > 0, 1, 0)
    }
    if (terms == "binary" || terms == "quadratic") {
        J <- tcrossprod(y, x)
        A <- rowSums(x^2)[col(J)]
        B <- rowSums(y^2)[row(J)]
    }
    ## terms = "minimum" should be written in C like in designdist
    if (terms == "minimum") {
        J <- matrix(0, Ny, Nx)
        for (j in seq_len(Nx))
            for (i in seq_len(Ny))
                J[i,j] <- sum(pmin.int(y[i,], x[j,]))
        dimnames(J) <- list(rownames(y), rownames(x))
        A <- rowSums(x)[col(J)]
        B <- rowSums(y)[row(J)]
    }
    Jattr <- attributes(J)
    J <- as.vector(J)
    ## 2x2 contingency table notation
    if (abcd) {
        a <- J
        b <- A - J
        c <- B - J
        d <- P - A - B + J
    }
    ## beta diversity notation
    if (alphagamma) {
        alpha <- (A + B)/2
        gamma <- A + B - J
        delta <- abs(A - B)/2
    }
    dis <- eval(parse(text = method))
    attributes(dis) <- Jattr
    attr(dis, "call") <- match.call()
    if (missing(name))
        attr(dis, "method") <- paste(terms, method)
    else attr(dis, "method") <- name
    if (!missing(maxdist)) {
        if (!is.na(maxdist) && any(dis > maxdist)) {
            warning("'maxdist' was lower than some distances: setting to NA")
            maxdist <- NA
        }
        attr(dis, "maxdist") <- maxdist
    }
    dis
}

## similar to designdist, but uses Chao's terms U & V instead of J, A,
## B (or their derived terms) in designdist. I considered having this
## as an option 'terms = "chao"' in designdist, but there really is so
## little in common and too many if's needed.

`chaodist` <-
    function(x, method = "1 - 2*U*V/(U+V)", name)
{
    x <- as.matrix(x)
    ## need integer data
    if (!isTRUE(all.equal(x, round(x))))
        stop("function accepts only integers (counts)")
    x <- round(x) # to be sure since as.integer(sqrt(3)^2) == 2
    N <- nrow(x)
    ## do_chaoterms returns a list with U, V which are non-classed
    ## vectors where the order of terms matches 'dist' objects
    vu <- .Call(do_chaoterms, x)
    U <- vu$U
    V <- vu$V
    ## dissimilarities
    dis <- eval(parse(text = method))
    dis <- structure(dis, Size = N, Labels = rownames(x), Diag = FALSE,
                     Upper = FALSE, call = match.call(), class = "dist")
    if (missing(name))
        attr(dis, "method") <- paste("chao", method)
    else
        attr(dis, "method") <- name
    dis
}

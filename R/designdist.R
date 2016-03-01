`designdist` <-
    function (x, method = "(A+B-2*J)/(A+B)",
              terms = c("binary", "quadratic", "minimum"),
              abcd = FALSE, alphagamma = FALSE, name) 
{
    terms <- match.arg(terms)
    if ((abcd || alphagamma) && terms != "binary")
        warning("Perhaps terms should be 'binary' with 'abcd' or 'alphagamma'?")
    x <- as.matrix(x)
    N <- nrow(x)
    P <- ncol(x)
    if (terms == "binary") 
        x <- ifelse(x > 0, 1, 0)
    if (terms == "binary" || terms == "quadratic") 
        x <- tcrossprod(x)
    if (terms == "minimum") {
        r <- rowSums(x)
        x <- dist(x, "manhattan")
        x <- (outer(r, r, "+") - as.matrix(x))/2
    }
    d <- diag(x)
    A <- as.dist(outer(rep(1, N), d))
    B <- as.dist(outer(d, rep(1, N)))
    J <- as.dist(x)
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
    dis
}

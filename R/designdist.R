`designdist` <-
    function (x, method = "(A+B-2*J)/(A+B)",
              terms = c("binary", "quadratic", "minimum"),
              abcd = FALSE, name) 
{
    terms <- match.arg(terms)
    if (abcd && terms != "binary")
        warning("abcd = TRUE and terms are not 'binary':\nresults may be meaningless")
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
    if (abcd) {
        a <- J
        b <- A - J
        c <- B - J
        d <- P - A - B + J
    }
    dis <- eval(parse(text = method))
    attributes(dis) <- attributes(J)
    attr(dis, "call") <- match.call()
    if (missing(name)) 
        attr(dis, "method") <- paste(terms, method)
    else attr(dis, "method") <- name
    dis
}

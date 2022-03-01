`betadiver` <-
    function(x, method = NA, order = FALSE, help = FALSE,  ...)
{
    beta <- list("w"="(b+c)/(2*a+b+c)", "-1"="(b+c)/(2*a+b+c)", "c"="(b+c)/2",
                 "wb"="b+c", "r"="2*b*c/((a+b+c)^2-2*b*c)",
                 "I"="log(2*a+b+c) - 2*a*log(2)/(2*a+b+c) - ((a+b)*log(a+b) + (a+c)*log(a+c)) / (2*a+b+c)",
                 "e"="exp(log(2*a+b+c) - 2*a*log(2)/(2*a+b+c) - ((a+b)*log(a+b) + (a+c)*log(a+c)) / (2*a+b+c))-1",
                 "t"="(b+c)/(2*a+b+c)", "me"="(b+c)/(2*a+b+c)",
                 "j"="a/(a+b+c)", "sor"="2*a/(2*a+b+c)",
                 "m"="(2*a+b+c)*(b+c)/(a+b+c)",
                 "-2"="pmin(b,c)/(pmax(b,c)+a)",
                 "co"="(a*c+a*b+2*b*c)/(2*(a+b)*(a+c))",
                 "cc"="(b+c)/(a+b+c)", "g"="(b+c)/(a+b+c)",
                 "-3"="pmin(b,c)/(a+b+c)", "l"="(b+c)/2",
                 "19"="2*(b*c+1)/(a+b+c)/(a+b+c-1)",
                 "hk"="(b+c)/(2*a+b+c)", "rlb"="a/(a+c)",
                 "sim"="pmin(b,c)/(pmin(b,c)+a)",
                 "gl"="2*abs(b-c)/(2*a+b+c)",
                 "z"="(log(2)-log(2*a+b+c)+log(a+b+c))/log(2)"
                 )
    if (help) {
        for (i in seq_along(beta))
            writeLines(strwrap(paste(i, " \"", names(beta[i]),
                                     "\" = ", beta[[i]], "\n", sep="")))
        return(invisible(NULL))
    }
    x <- ifelse(x > 0, 1, 0)
    if (order) {
        x <- x[order(rowSums(x)),]
    }
    d <- tcrossprod(x)
    a <- as.dist(d)
    S <- diag(d)
    N <- length(S)
    b <- as.dist(matrix(rep(S, N), nrow=N)) - a
    c <- as.dist(matrix(rep(S, each=N), nrow=N)) - a
    if (is.na(method) || is.null(method) || is.logical(method) && !method) {
        out <- list(a = a, b = b, c = c)
        class(out) <- "betadiver"
        return(out)
    }
    out <- eval(parse(text=beta[[method]]))
    out <- as.dist(out)
    attr(out, "method") <- paste("beta", names(beta[method]), sep=".")
    attr(out, "call") <- match.call()
    out
}

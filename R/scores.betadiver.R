`scores.betadiver` <-
    function(x, triangular = TRUE,  ...)
{
    if (triangular) {
        tot <- x$a + x$b + x$c
        a <- x$a/tot
        c <- x$c/tot
        y <- sqrt(0.75)*a
        x <- c + a/2
        out <- cbind(x, y)
    } else {
        out <- sapply(x, cbind)
    }
    out
}

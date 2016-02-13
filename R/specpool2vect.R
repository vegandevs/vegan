"specpool2vect" <-
function(X, index = c("jack1","jack2", "chao", "boot", "Species"))
{
    pool <- attr(X, "pool")
    index <- match.arg(index)
    X[[index]][pool]
}

"specpool2vect" <-
function(X, index = c("Jack.1","Jack.2", "Chao", "Boot", "Species"))
{
    pool <- attr(X, "pool")
    index <- match.arg(index)
    X[[index]][pool]
}

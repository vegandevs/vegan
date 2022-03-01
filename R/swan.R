swan <-
function (x, maxit = Inf, type = 0)
{
    zeros <- -Inf
    iter <- 0
    while(zeros != (zeros <- sum(x == 0)) && any(x == 0) &&
          iter < maxit) {
        x[x > 0] <- x[x > 0] - min(x[x > 0]) + 1
        x[x == 0] <- beals(x, type = type)[x == 0]
        iter <- iter + 1
    }
    x
}
### (Ecology 51, 89-102; 1970).

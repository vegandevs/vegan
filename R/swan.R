swan <-
function (x)
{
    while(any(x == 0)) {
        x[x > 0] <- x[x > 0] - min(x[x > 0]) + 1
        x[x == 0] <- beals(x)[x == 0]
    }
    x
}
### (Ecology 51, 89-102; 1970).

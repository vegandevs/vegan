rrarefy <-
function(x, sample)
{
    if (length(sample) > 1 && length(sample) != nrow(x))
        stop("length of 'sample' and number of rows of 'x' do not match")
    sample <- rep(sample, length=nrow(x))
    colnames(x) <- colnames(x, do.NULL = FALSE)
    nm <- colnames(x)
    for (i in 1:nrow(x)) {
        row <- sample(rep(nm, times=x[i,]), sample[i])
        row <- table(row)
        ind <- names(row)
        x[i,] <- 0
        x[i,ind] <- row
    }
    x
}

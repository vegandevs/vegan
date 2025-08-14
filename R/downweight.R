### Support functions for decorana

## Hill's downweighting

## An exported function that can be called outside decorana

`downweight` <-
    function (veg, fraction = 5)
{
    Const1 <- 1e-10
    if (fraction < 1)
        fraction <- 1/fraction
    veg <- as.matrix(veg, rownames.force = TRUE)
    yeig1 <- colSums(veg)
    y2 <- colSums(veg^2) + Const1
    y2 <- yeig1^2/y2
    amax <- max(y2)/fraction
    v <- rep(1, ncol(veg))
    downers <- y2 < amax
    v[downers] <- (y2/amax)[downers]
    veg <- sweep(veg, 2, v, "*")
    attr(veg, "v") <- v
    attr(veg, "fraction") <- fraction
    veg
}

## Hill's piecewise tranformation. Values of before are replaced with
## values of after, and intermediary values with linear interpolation.

## Not exported: if you think you need something like this, find a
## better tool in R.

`beforeafter` <-
    function(x, before, after)
{
    if (is.null(before) || is.null(after))
        stop("both 'before' and 'after' must be given", call. = FALSE)
    if (is.unsorted(before))
        stop("'before' must be sorted", call. = FALSE)
    if (length(before) != length(after))
        stop("'before' and 'after' must have same lengths", call. = FALSE)
    for(i in seq_len(nrow(x))) {
        k <- x[i,] > 0
        x[i, k] <- approx(before, after, x[i, k], rule = 2)$y
    }
    x
}

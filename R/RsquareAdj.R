"RsquareAdj" <-
function(x, n, m, ...)
{
    if (m >= (n-1))
        NA
    else
        1 - (1-x)*(n-1)/(n-m-1)
}


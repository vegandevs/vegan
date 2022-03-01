### Internal function for double centring of a *matrix* of
### dissimilarities. We used .C("dblcen", ..., PACKAGE = "stats")
### which does not dublicate its argument, but it was removed from R
### in r60360 | ripley | 2012-08-22 07:59:00 UTC (Wed, 22 Aug 2012)
### "more conversion to .Call, clean up". Input 'x' *must* be a
### matrix. This was originally an internal function in betadisper.R
### (commit 7cbd4529 Thu Aug 23 08:45:31 2012 +0000)
GowerDblcen <- function(x, na.rm = TRUE)
{
    cnt <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2L, cnt, check.margin = FALSE)
    cnt <- rowMeans(x, na.rm = na.rm)
    sweep(x, 1L, cnt, check.margin = FALSE)
}

### Internal functions to find additive constants to non-diagonal
### dissimilarities so that there are no negative eigenvalues. The
### Cailliez constant is added to dissimilarities and the Lingoes
### constant is added to squared dissimilarities. Legendre & Anderson
### (Ecol Monogr 69, 1-24; 1999) recommend Lingoes, but
### stats::cmdscale() only provides Cailliez. Input parameters: d are
### a matrix of dissimilarities.

addCailliez <- function(d)
{
    n <- nrow(d)
    q1 <- seq_len(n)
    q2 <- n + q1
    ## Cailliez makes a 2x2 block matrix with blocks of n x n elements.
    ## Blocks anti-clockwise, upper left [0]
    z <- matrix(0, 2*n, 2*n)
    diag(z[q2,q1]) <- -1
    z[q1,q2] <- -GowerDblcen(d^2)
    z[q2,q2] <- GowerDblcen(2 * d)
    ## Largest real eigenvalue
    e <- eigen(z, symmetric = FALSE, only.values = TRUE)$values
    out <- max(Re(e))
    max(out, 0)
}

addLingoes <- function(d)
{
    ## smallest negative eigenvalue (or zero)
    d <- -GowerDblcen(d^2)/2
    e <- eigen(d, symmetric = TRUE, only.values = TRUE)$values
    out <- min(e)
    max(-out, 0)
}

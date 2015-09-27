### Internal function for double centring of a *matrix* of
### dissimilarities. We used .C("dblcen", ..., PACKAGE = "stats")
### which does not dublicate its argument, but it was removed from R
### in r60360 | ripley | 2012-08-22 07:59:00 UTC (Wed, 22 Aug 2012)
### "more conversion to .Call, clean up". Input 'x' *must* be a
### matrix. This was originally an internal function in betadisper.R
### (commit 7cbd4529 Thu Aug 23 08:45:31 2012 +0000)
GowerDblcen <- function(x, na.rm = TRUE) {
    cnt <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2L, cnt, check.margin = FALSE)
    cnt <- rowMeans(x, na.rm = na.rm)
    sweep(x, 1L, cnt, check.margin = FALSE)
}

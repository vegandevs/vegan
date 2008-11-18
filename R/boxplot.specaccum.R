"boxplot.specaccum" <-
    function(x, add=FALSE, ...)
{
    if (x$method != "random")
        stop("boxplot available only for method=\"random\"")
    if (!add) {
        plot(x$sites, x$richness, type="n", xlab="Sites", ylab="Species",
             ylim=c(1, max(x$richness)),  ...)
    }
    tmp <- boxplot(data.frame(t(x$perm)), add=TRUE, at=x$sites, axes=FALSE, ...)
    invisible(tmp)
}

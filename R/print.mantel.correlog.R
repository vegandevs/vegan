'print.mantel.correlog' <- function(x, ...)
{
    cat('\nMantel Correlogram Analysis\n')
    cat('\nCall:\n','\n')
    cat(deparse(x$call),'\n')
    cat('\n')
	printCoefmat(x$mantel.res, P.values=TRUE, signif.stars=TRUE, Pvalues = TRUE)
    invisible(x) 
}
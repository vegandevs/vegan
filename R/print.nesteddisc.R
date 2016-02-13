"print.nesteddisc" <-
function(x, ...)
{
    cat("nestedness discrepancy:", x$statistic, "\n")
    if(x$ties)
        cat("There are tied column frequencies: result can depend on input order\n")
    invisible(x)
}

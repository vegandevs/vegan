`print.summary.meandist` <-
    function(x, ...)
{
    cat("\nMean distances:\n")
    tab <- rbind("within groups" = x$W,
                 "between groups" = x$B,
                 "overall" = x$D)
    colnames(tab) <- "Average"
    print(tab, ...)
    cat("\nSummary statistics:\n")
    tab <- rbind("MRPP A weights n" = x$A1,
                 "MRPP A weights n-1" = x$A2,
                 "MRPP A weights n(n-1)"= x$A3,
                 "Classification strength"=x$CS)
    colnames(tab) <- "Statistic"
    print(tab, ...)
    invisible(x)
}


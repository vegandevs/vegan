`print.taxondive` <-
    function (x, ...) 
{
    out <- cbind(x$Species, x$D, x$Dstar, x$Lambda, x$Dplus, x$SDplus)
    out <- rbind(out, Expected = c(NA, x$ED, x$EDstar, NA, x$EDplus, NA))
    colnames(out) <- c("Species", "Delta", "Delta*", "Lambda+", 
                       "Delta+", "S Delta+")
    printCoefmat(out, na.print = "")
    invisible(x)
}

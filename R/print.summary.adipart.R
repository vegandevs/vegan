print.summary.adipart <-
function(x, ...)
{
    cat("Additive diversity partitioning summary\n\nCall: ")
    print(x$call)
    cat("\n")
    print(x$divres, quote=FALSE)
    cat("---\nSignif. codes:  0 \"***\" 0.001 \"**\" 0.01 \"*\" 0.05 \".\" 0.1 \"ns\" 1\n")
    cat("Based on", x$times, "permutations.\n")
}

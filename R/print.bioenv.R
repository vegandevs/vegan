"print.bioenv" <-
    function (x, ...) 
{
    cat("\nCall:\n")
    cat(deparse(x$call), "\n")
    cat("\nSubset of environmental variables with best correlation to community data.\n\n")
    cat("Correlations:     ", x$method, "\n")
    cat("Dissimilarities:  ", x$index, "\n\n")
    i <- which.max(lapply(x$models, function(tmp) tmp$est))
    cat("Best model has", i, "parameters (max.", x$upto, "allowed):\n")
    cat(paste(x$names[x$models[[i]]$best], collapse = " "))
    cat("\nwith correlation ", x$models[[i]]$est, "\n\n")
    invisible(x)
}

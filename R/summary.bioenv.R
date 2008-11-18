"summary.bioenv" <-
    function(object, ...)
{
    x <- object$models
    nam <- object$names
    size <- seq(1:length(x))
    cor <- unlist(lapply(x, function(tmp) tmp$est))
    pars <- unlist(lapply(x, function(tmp) paste(nam[tmp$best], collapse=" ")))
    out <- list(size = size, correlation = cor, variables = pars)
    class(out) <- "summary.bioenv"
    out
}

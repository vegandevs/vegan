## this is function to create a commsim object, does some checks
## there is a finite number of useful arguments here
## but I added ... to allow for unforeseen algorithms,
## or being able to reference to external objects
commsim <- 
function(method, fun, binary, isSeq, mode) 
{
    fun <- if (!missing(fun))
        match.fun(fun) else stop("'fun' missing")
    if (any(!(names(formals(fun)) %in% 
        c("x", "n", "nr", "nc", "rs", "cs", "rf", "cf", "s", "fill", "thin", "..."))))
            stop("unexpected arguments in 'fun'")
    out <- structure(list(method = if (!missing(method))
            as.character(method)[1L] else stop("'method' missing"),
        binary = if (!missing(binary))
            as.logical(binary)[1L] else stop("'binary' missing"),
        isSeq = if (!missing(isSeq))
            as.logical(isSeq)[1L] else stop("'isSeq' missing"),
        mode = if (!missing(mode))
            match.arg(as.character(mode)[1L],
            c("integer", "double")) else stop("'mode' missing"),
        fun = fun), class = "commsim")
    out
}

`print.allPerms` <- function(x, ...) {
    dims <- dim(x)
    control <- attr(x, "control")
    observed <- attr(x, "observed")
    attributes(x) <- NULL
    dim(x) <- dims
    print(x)
    return(invisible(x))
}

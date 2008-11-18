`print.permControl` <- function(x, ...)
{
    ## only for objects of correct class
    stopifnot(class(x) == "permControl")
    ## set-up the messages we wish to print
    if (!is.null(x$strata)) {
        if(x$permute.strata)
            msg.strata <- paste("Permutations stratified between '",
                                x$name.strata, "'\n", sep = "")
        else
            msg.strata <- paste("Permutations stratified within '",
                                x$name.strata, "'\n", sep = "")
    } else {
        msg.strata <- "Permutations are unstratified\n"
    }
    msg.ptype <- paste("Permutation type:", x$type, "\n")
    if(x$type == "grid")
        msg.grid <- paste("Data are spatial grid(s) of dimension",
                          x$nrow, "*", x$ncol, "\n")
    msg.nperm <- paste("No. of permutations:", x$nperm,
                       ifelse(x$complete, "(complete enumeration)", ""),
                       "\n")
    msg.mirror <- paste("Mirrored permutations?:",
                        ifelse(x$mirror, "Yes", "No"), "\n")
    msg.constant <- paste("Use same permutation within strata?:",
                          ifelse(x$constant, "Yes", "No"), "\n")
    ## print out the messages
    cat("\n")
    cat(msg.nperm)
    cat(msg.ptype)
    if(exists("msg.grid"))
        cat(msg.grid)
    cat(msg.strata)
    cat(msg.mirror)
    cat(msg.constant)
    cat("\n")
}

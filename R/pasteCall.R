`pasteCall` <- function (call, prefix = "Call:") 
{
    call.str <- paste(deparse(call), collapse = " ")
    call.str <- paste(prefix, call.str, "\n", sep = " ")
    return(call.str)
}

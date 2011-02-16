`pasteCall` <- function (call, prefix = "Call:") 
{
    call.str <- paste(deparse(call), collapse = " ")
    paste(prefix, call.str, "\n", sep = " ")
}

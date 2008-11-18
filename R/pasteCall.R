`pasteCall` <- function (call, prefix = "Call:") 
{
    objcall <- deparse(call)
    if (length(objcall) > 1) {
        call.str <- ""
        for (i in seq(along = objcall)) {
            call.str <- paste(call.str, objcall[i], sep = " ")
        }
    }
    else {
        call.str <- objcall
    }
    call.str <- paste(prefix, call.str, "\n", sep = " ")
    return(call.str)
}

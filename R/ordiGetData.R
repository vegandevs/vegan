`ordiGetData` <-
function (call, env) 
{
    call$scale <- call$distance <- call$comm <- call$add <- NULL
    call$na.action <- na.pass
    call[[2]] <- NULL
    call[[1]] <- as.name("model.frame")
    eval(call, env)
}


`ordiGetData` <-
function (call, env) 
{
    call$scale <- call$distance <- call$comm <- call$add <-
        call$dfun <- call$sqrt.dist <- call$metaMDSdist <- call$subset <- NULL
    call$na.action <- na.pass
    call[[2]] <- NULL
    call[[1]] <- as.name("model.frame")
    eval(call, env)
}

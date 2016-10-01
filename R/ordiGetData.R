`ordiGetData` <-
function (call, env)
{
    eval(call$data, env, enclos = .GlobalEnv)
}

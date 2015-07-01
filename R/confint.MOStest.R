`confint.MOStest` <-
    function (object, parm = 1, level = 0.95, ...) 
{
    confint(profile(object), level = level, ...)
}

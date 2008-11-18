`confint.MOStest` <-
    function (object, parm = 1, level = 0.95, ...) 
{
    require(MASS) || stop("requires packages MASS")
    confint(profile(object), level = level, ...)
}

"confint.fisherfit" <-
    function (object, parm, level=0.95,  ...)
{
    if (!require(MASS)) stop("Needs packages MASS .. not found")
    confint(profile(object), level=level, ...)
}

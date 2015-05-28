### density and densityplot were deprecated in vegan 2.2-0 in favour
### of corresponding functions for permustats()

`density.anosim` <-
    function(x, ...)
{
    .Defunct("densityplot(permustats(<anosim.result>))",
                package="vegan")
}

`density.adonis` <-
    function(x, ...)
{
    .Defunct("densityplot(permustats(<adonis.result>))",
                package="vegan")
}

`densityplot.adonis` <-
    function(x, data, xlab = "Null", ...)
{
    .Defunct("densityplot(permustats(<adonis.result>))",
                package="vegan")
}

`density.mantel` <-
    function(x, ...)
{
    .Defunct("densityplot(permustats(<mantel.result>))",
                package="vegan")
}

`density.mrpp` <-
    function(x, ...)
{
    .Defunct("densityplot(permustats(<mrpp.result>))",
                package="vegan")
}

`density.permutest.cca` <-
    function(x, ...)
{
    .Defunct("densityplot(permustats(<permutest.result>))",
                package="vegan")
}

`density.protest` <-
    function(x, ...)
{
    .Defunct("densityplot(permustats(<protest.result>))",
                package="vegan")
}

`plot.vegandensity` <-
    function (x, main = NULL, xlab = NULL, ylab = "Density", type = "l", 
    zero.line = TRUE, obs.line = TRUE, ...) 
{
    .Defunct("permustats methods", package = "vegan")
}

`density.oecosimu` <-
    function(x, ...)
{
    .Defunct("densityplot(permustats(<oecosimu.result>))",
                package="vegan") 
}

`densityplot.oecosimu` <-
    function(x, data, xlab = "Simulated", ...)
{
    .Defunct("densityplot(permustats(<oecosimu.result>))",
                package="vegan")
}

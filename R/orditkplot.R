###
### Editable Tcl/Tk plot for ordination
###
`orditkplot` <-
    function(...)
{
    if (requireNamespace("vegan3d")) {
        message("orditkplot was moved to vegan3d and is deprecated in vegan")
        vegan3d::orditkplot(...)
    } else {
        stop("orditkplot was moved to CRAN package vegan3d:
              install vegan3d from CRAN and use vegan3d::orditkplot")
    }
}

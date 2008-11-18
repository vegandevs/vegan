`ordicloud` <-
    function(x, data = NULL, formula, display = "sites", choices=1:3,
             panel = "panel.ordi3d",
             prepanel = "prepanel.ordi3d", ...)
{
  localCloud <- function(..., shrink, origin, scaling) cloud(...)
  require(lattice) || stop("requires package 'lattice'")
  x <- as.data.frame(scores(x, display = display, choices = choices, ...))
  if (!is.null(data))
    x <- cbind(x, data)
  if (missing(formula)) {
    v <- colnames(x)
    formula <- as.formula(paste(v[2], "~", v[1], "*", v[3]))
  }
  localCloud(formula, data = x, panel = panel,  prepanel = prepanel, ...)
}

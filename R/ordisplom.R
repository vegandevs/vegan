`ordisplom` <-
    function(x, data=NULL, formula = NULL,  display = "sites", choices = 1:3,
             panel = "panel.ordi", type = "p", ...)
{
  localSplom <- function(..., shrink, origin, scaling) splom(...)
  x <- as.data.frame(scores(x, display = display, choices = choices, ...))
  if (is.null(data))
    data <- x
  else if (is.null(formula))
    x <- cbind(x, data)
  ## type = "biplot" is not (yet?) implemented
  env <- list(arrows = NULL, centres = NULL)
  if (is.null(formula))
    pl <- localSplom(x, panel = panel, type = type, biplot = env, ...)
  else {
    formula <- as.formula(gsub("\\.", "x", deparse(formula)))
    pl <- localSplom(x = formula, data = data,  panel = panel, type = type,
                     bitplot = env, ...)
  }
  pl
}

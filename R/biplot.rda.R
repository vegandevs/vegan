## biplot.rda
##
## draws pca biplots with species as arrows
##

`biplot.cca` <-
    function(x, ...)
{
    if (!inherits(x, "rda"))
        stop("biplot can be used only with linear ordination (e.g., PCA)")
    else
        NextMethod("biplot")
}

`biplot.rda` <- function(x, choices = c(1, 2), scaling = "species",
                         display = c("sites", "species"),
                         type, xlim, ylim, col = c(1,2), const,
                         correlation = FALSE, ...) {
  if(!inherits(x, "rda"))
      stop("'biplot.rda' is only for objects of class 'rda'")
  if(!is.null(x$CCA))
      stop("'biplot.rda' not suitable for models with constraints")
  TYPES <- c("text", "points", "none")
  display <- match.arg(display, several.ok = TRUE)
  if (length(col) == 1)
      col <- c(col,col)
  spe.par <- list(col = col[2], arrows = TRUE)
  sit.par <- list(col = col[1])
  if (!missing(type)) {
      type <- match.arg(type, TYPES, several.ok = TRUE)
      if(length(type) < 2)
          type <- rep(type, 2)
      spe.par <- modifyList(spe.par, list(type = type[1]))
      sit.par <- modifyList(sit.par, list(type = type[2]))
  }
  plot(x, choices = choices, scaling = scaling, display = display,
       xlim = xlim, ylim = ylim, const = const, correlation = correlation,
       spe.par = spe.par, sit.par = sit.par, ...)
}

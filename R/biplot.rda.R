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
        NextMethod("biplot", x, ...)
}

`biplot.rda` <- function(x, choices = c(1, 2), scaling = 2,
                       display = c("sites", "species"),
                       type, xlim, ylim, col = c(1,2), const, ...) {
  if(!inherits(x, "rda"))
      stop("'biplot.rda' is only for objects of class 'rda'")
  if(!is.null(x$CCA))
      stop("'biplot.rda' not suitable for models with constraints")
  TYPES <- c("text", "points", "none")
  display <- match.arg(display, several.ok = TRUE)
  if (length(col) == 1)
      col <- c(col,col)
  g <- scores(x, choices = choices, display = display,
              scaling = scaling, const)
  if (!is.list(g)) {
      g <- list(default = g)
      names(g) <- display
  }
  if (missing(type)) {
      nitlimit <- 80
      nit <- max(nrow(g$species), nrow(g$sites))
      if (nit > nitlimit)
          type <- rep("points", 2)
      else type <- rep("text", 2)
  }
  else type <- match.arg(type, TYPES, several.ok = TRUE)
  if(length(type) < 2)
      type <- rep(type, 2)
  if (missing(xlim))
      xlim <- range(g$species[, 1], g$sites[, 1], na.rm = TRUE)
  if (missing(ylim))
      ylim <- range(g$species[, 2], g$sites[, 2], na.rm = TRUE)
  plot(g[[1]], xlim = xlim, ylim = ylim, type = "n", asp = 1,
       ...)
  abline(h = 0, lty = 3)
  abline(v = 0, lty = 3)
  if (!is.null(g$species)) {
      if (type[1] == "points")
          arrlen <- 1
      else
          arrlen <- 0.85
      if (type[1] != "none")
          arrows(0, 0, g$species[,1] * arrlen, g$species[, 2] * arrlen,
                 col = col[2], length = 0.05)
      if (type[1] == "text")
          text(g$species, rownames(g$species),
               col = col[2], cex = 0.7)
  }
  if (!is.null(g$sites)) {
      if (type[2] == "text")
          text(g$sites, rownames(g$sites), cex = 0.7, col = col[1])
      else if (type[2] == "points")
          points(g$sites, pch = 1, cex = 0.7, col = col[1])
  }
  class(g) <- "ordiplot"
  invisible(g)
}

`ordixyplot` <-
    function(x, data = NULL, formula, display = "sites", choices=1:3,
             panel = "panel.ordi", aspect = "iso", envfit,
             type = c("p", "biplot"),  ...)
{
  localXyplot <- function(..., shrink, origin, scaling) xyplot(...)
  require(lattice) || stop("requires package 'lattice'")
  p <- as.data.frame(scores(x, display = display, choices = choices, ...))
  if (!is.null(data))
    p <- cbind(p, data)
    if (missing(formula)) {
      v <- colnames(p)
      formula <- as.formula(paste(v[2], "~", v[1]))
    }
  if ("biplot" %in% type && ((!is.null(x$CCA) && x$CCA$rank > 0) ||
                             !missing(envfit))) {
    if (missing(envfit))
      envfit <- NULL
    env <- ordilattice.getEnvfit(formula, x, envfit, choices, ...)
    if (!is.null(env$arrows)) {
      mul <- apply(p[,colnames(env$arrows)], 2, range)/apply(env$arrows, 2, range)
      mul <- min(mul[is.finite(mul) & mul > 0])
            env$arrows <- mul * env$arrows
    }
  } else {
    env <- NULL
  }
  localXyplot(formula, data = p, panel = panel,  aspect = aspect, biplot = env, type = type,
         ...)
}

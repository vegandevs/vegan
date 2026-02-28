\name{vegan-deprecated}

\alias{vegan-deprecated}
\alias{as.mcmc.oecosimu}
\alias{as.mcmc.permat}
\alias{ordixyplot}
\alias{panel.ordi}
\alias{panel.ordiarrows}
\alias{ordilattice.getEnvfit}
\alias{plot.poolaccum}
\alias{plot.renyiaccum}
\alias{plot.renyi}
\alias{densityplot.permustats}
\alias{qqmath.permustats}
\alias{permulattice}

%------ NOTE:  ../R/vegan-deprecated.R   must be synchronized with this!
\title{Deprecated Functions in vegan package}
%------ PLEASE: one \alias{.} for EACH ! (+ one \usage{} & \arguments{} for all)

\description{
  These functions are provided for compatibility with older versions of
  \pkg{vegan} only, and may be defunct as soon as the next release.
}

\usage{
## Following lattice functions are deprecated. Use ggplot2 functions
## in the CRAN package ggvegan instead

ordixyplot(x, data = NULL, formula, display = "sites", choices = 1:3,
    panel = "panel.ordi", aspect = "iso", envfit,
    type = c("p", "biplot"), ...)
\method{plot}{poolaccum}(x, alpha = 0.05, type = c("l","g"), ...)
\method{plot}{renyi}(x, ...)
\method{plot}{renyiaccum}(x, what = c("Collector", "mean", "Qnt 0.025",
    "Qnt 0.975"),
     type = "l", ...)
permulattice(x, plot = c("densityplot", "qqmath"), observed = TRUE,
    axislab = "Permutations", ...)
\method{densityplot}{permustats}(x, data, observed = TRUE,
    xlab = "Permutations", ...)
\method{qqmath}{permustats}(x, data, observed = TRUE, sd.scale = FALSE,
    ylab = "Permutations", ...)

## use toCoda instead
as.mcmc.oecosimu(x)
as.mcmc.permat(x)
}

\arguments{

  \item{x}{input object}
  \item{data}{ Optional data to amend ordination results. The ordination
    results are found from \code{x}, but you may give here data for other
    variables needed in plots. Typically these are environmental data.}
  \item{formula}{ Formula to define the plots. A default formula will be
    used if this is omitted. The
    ordination axes must be called by the same names as in the
    ordination results (and these names vary among methods). }
  \item{display}{ The kind of scores: an argument passed to
    \code{\link{scores}}. }
  \item{choices}{ The axes selected: an argument passed to
    \code{\link{scores}}. } 
  \item{panel}{ The name of the panel function. } 
  \item{aspect}{The aspect of the plot (passed to the \pkg{lattice}
    function).}
  \item{envfit}{Result of \code{\link{envfit}} function displayed in
    \code{ordixyplot}. Please note that this needs same \code{choices}
    as \code{ordixyplot}.}
  \item{type}{The type of plot. This knows the same alternatives as
    \code{\link[lattice]{panel.xyplot}}. In addition \code{ordixyplot}
    has alternatives \code{"biplot"}, \code{"arrows"} and
    \code{"polygon"}. The first displays fitted vectors and factor
    centroids of \code{envfit}, or in constrained ordination, the biplot
    arrows and factor centroids if \code{envfit} is not given. The
    second (\code{type = "arrows"}) is a trellis variant of
    \code{\link{ordiarrows}} and draws arrows by \code{groups}. The line
    parameters are controlled by \code{\link[lattice]{trellis.par.set}}
    for \code{superpose.line}, and the user can set \code{length},
    \code{angle} and \code{ends} parameters of
    \code{\link[lattice]{panel.arrows}}.  The last one (\code{type =
    "polygon"}) draws a polygon enclosing all points in a panel over a
    polygon enclosing all points in the data. The overall polygon is
    controlled by Trellis parameters \code{\link[lattice]{trellis.par.set}}
    \code{plot.polygon} and \code{superpose.polygon}.}
  \item{what}{Items to be plotted.}
    \item{plot}{Use \pkg{lattice} function
    \code{\link[lattice]{densityplot}} or \code{\link[lattice]{qqmath}}.}
  \item{xlab, ylab, axislab}{Label for the axis displaying permutation
    values.}
  \item{observed}{Add observed statistic among permutations.}
  \item{sd.scale}{Scale permutations to unit standard deviation and observed
    statistic to standardized effect size.}
  \item{alpha}{Level of quantiles shown. This proportion will be left outside
    symmetric limits.}
  \item{\dots}{ Arguments passed to \code{\link{scores}} methods or
    \pkg{lattice} functions.}
}

\details{

  \code{as.mcmc} functions were replaced with \code{\link{toCoda}}.

  Better alternatives to \pkg{lattice} functions are provided by
  \pkg{ggplot2} functions in the \CRANpkg{ggvegan}
  package. \code{ggvegan::autoplot} functions are most similar to these
  deprecated functions.
}

\seealso{
  \code{\link{Deprecated}}
}
\keyword{misc}

\name{vegan-deprecated-lattice}
\alias{ordicloud}
\alias{ordisplom}
\alias{panel.ordi3d}
\alias{prepanel.ordi3d}
\alias{ordiresids}
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

\title{Deprecated lattice Functions in vegan}

\description{Lattice functions are going to be deprecated and removed
  from \pkg{vegan}. They are replaced with \pkg{ggplot2} functions in
  CRAN package \CRANpkg{ggvegan}.}

\usage{
ordicloud(x, data = NULL, formula, display = "sites", choices = 1:3,
    panel = "panel.ordi3d", prepanel = "prepanel.ordi3d", ...)
ordisplom(x, data = NULL, formula = NULL, display = "sites", choices = 1:3, 
    panel = "panel.ordi", type = "p", ...)
ordiresids(x, kind = c("residuals", "scale", "qqmath"),
    residuals = "working", type = c("p", "smooth", "g"),
    formula, ...)
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

}

\arguments{
  \item{x}{Input object.}
  \item{kind}{The type of plot: residuals or absolute values of
    residuals against fitted values, or quantile plot of residuals with
    \code{\link[lattice]{qqmath}}.}
  \item{residuals}{The type of residuals with choices \code{"working"},
    \code{"response"}, \code{"standardized"} and \code{"studentized"}.}
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
  \item{panel, prepanel}{ The name of the panel or prepanel function. } 
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

  Trellis (or \pkg{lattice}) functions were added to \pkg{vegan}
  mostly in 2008 to 2009. In that time they were the only alternative
  of the kind, but now there are better, more versatile and more
  user-friendly alternatives, mainly in \pkg{ggplot2}. CRAN package
  \pkg{ggvegan} provides modern alternatives to most lattice functions
  in \pkg{vegan}. The lattice functions in \pkg{vegan} will be
  deprecated as soon \pkg{ggvegan} provides a \pkg{ggplot2}
  alternative. The deprecated functions will be defunct in the next
  major release of \pkg{vegan}.

  The following functions are currently deprecated:

  \itemize{
    \item \code{ordicloud} was transferred to \CRANpkg{vegan3d} as
      \code{ordilattice3d}.
    \item \code{ordisplom} design is bad and deficient. If you want to
      have something similar, write your own code.
    \item \code{ordixyplot}: use \code{autoplot} or \code{ordiggplot}
      functions in \pkg{ggvegan}.
    \item \code{plot} functions for \code{\link{poolaccum}},
      \code{\link{renyi}} and \code{\link{renyiaccum}}: use \code{autoplot}
      in \pkg{ggvegan}.
    \item \code{permulattice}: use \code{autoplot} in \pkg{ggvegan}.
    \item \code{ordiresids} is not very useful, but you can directly
      access ordination results with \code{\link{fitted.cca}},
      \code{\link{residuals.cca}}, \code{\link{rstandard.cca}},
      \code{\link{rstudent.cca}} and other functions that were not
      available in \code{ordiresids}.
    }
}

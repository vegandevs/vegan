\name{vegan-deprecated-lattice}
\alias{ordicloud}
\alias{ordisplom}
\alias{panel.ordi3d}
\alias{prepanel.ordi3d}
\alias{ordiresids}

\title{Deprecated lattice Functions in vegan}

\description{Lattice functions are going to be deprecated and removed
  from vegan as soon as viable alternatives are available (or if they
  are unsatisfactory originally).}

\usage{
  ordicloud(x, data = NULL, formula, display = "sites", choices = 1:3,
    panel = "panel.ordi3d", prepanel = "prepanel.ordi3d", ...)
  ordisplom(x, data = NULL, formula = NULL, display = "sites", choices = 1:3, 
    panel = "panel.ordi", type = "p", ...)
  ordiresids(x, kind = c("residuals", "scale", "qqmath"),
    residuals = "working", type = c("p", "smooth", "g"),
    formula, ...)
}

\arguments{
  \item{x, data, formula, display, choices, panel, prepanel, type,
    \dots}{Similar parameters as in \code{\link{ordixyplot}}.}
  \item{kind}{The type of plot: residuals or absolute values of
    residuals against fitted values, or quantile plot of residuals with
    \code{\link[lattice]{qqmath}}.}
  \item{residuals}{The type of residuals with choices \code{"working"},
    \code{"response"}, \code{"standardized"} and \code{"studentized"}.}
}

\details{

  Trellis (or \CRANpkg{lattice}) functions were added to \pkg{vegan}
  mostly in 2008 to 2009. In that time they were the only alternative of
  the kind, but now there are better, more versatile and more
  user-friendly alternatives, mainly in \CRANpkg{ggplot2}. The lattice
  functions will be removed from \pkg{vegan} as soon as we can propose
  better modern alternatives.

  The following functions are currently deprecated:

  \itemize{
    \item \code{ordicloud} was transferred to \CRANpkg{vegan3d} as
      \code{ordilattice3d}.
    \item \code{ordisplom} design is bad and deficient. If you want to
      have something similar, write your own code.
    \item \code{ordiresids} is not very useful, but you can directly
      access ordination results with \code{\link{fitted.cca}},
      \code{\link{residuals.cca}}, \code{\link{rstandard.cca}},
      \code{\link{rstudent.cca}} and other functions that were not
      available in \code{ordiresids}.
    }
}

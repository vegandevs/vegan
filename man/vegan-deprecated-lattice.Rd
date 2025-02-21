\name{vegan-deprecated-lattice}
\alias{ordicloud}
\alias{ordisplom}
\alias{panel.ordi3d}
\alias{prepanel.ordi3d}

\title{Deprecated lattice Functions}

\description{Lattice functions are going to be deprecated and removed
  from vegan as soon as viable alternatives are available (or if they
  are unsatisfactory originally).}

\usage{
  ordicloud(x, data = NULL, formula, display = "sites", choices = 1:3,
    panel = "panel.ordi3d", prepanel = "prepanel.ordi3d", ...)
  ordisplom(x, data = NULL, formula = NULL, display = "sites", choices = 1:3, 
    panel = "panel.ordi", type = "p", ...)
}

\arguments{
  \item{x, data, formula, display, choices, panel, prepanel, type,
    \dots}{Similar parameters as in \code{\link{ordixyplot}}.}
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
    }
}
  
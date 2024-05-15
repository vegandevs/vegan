\name{vegan-deprecated}

\alias{orditkplot}
\alias{plot.orditkplot}
\alias{scores.orditkplot}
\alias{points.orditkplot}
\alias{text.orditkplot}
\alias{as.mcmc.oecosimu}
\alias{as.mcmc.permat}
\alias{vegan-deprecated}

%------ NOTE:  ../R/vegan-deprecated.R   must be synchronized with this!
\title{Deprecated Functions in vegan package}
%------ PLEASE: one \alias{.} for EACH ! (+ one \usage{} & \arguments{} for all)

\description{
  These functions are provided for compatibility with older versions of
  \pkg{vegan} only, and may be defunct as soon as the next release.
}

\usage{

## moved to vegan3d package: install from CRAN
orditkplot(x, display = "species", choices = 1:2, width, xlim, ylim, 
   tcex = 0.8, tcol, pch = 1,  pcol, pbg, pcex = 0.7, labels,  ...)
\method{plot}{orditkplot}(x, ...)
\method{points}{orditkplot}(x, pch = x$args$pch, cex = x$args$pcex,
       col = x$args$pcol, bg = x$args$pbg, ...)
\method{text}{orditkplot}(x, cex = x$args$tcex, col = x$args$tcol,
     font = attr(x$labels, "font"), ...)
\method{scores}{orditkplot}(x, display, ...)
## use toCoda instead
as.mcmc.oecosimu(x)
as.mcmc.permat(x)
}

\arguments{

  \item{x}{ordination result \code{orditkplot} or object to be
     tranformed \code{as.mcmc}.}
  \item{display, choices}{kind of scores and number of axes.}
  \item{width}{size of Tcl/Tk window in inches.}

  \item{xlim, ylim, tcex, tcol, pch, pcol, pbg, pcex, labels, cex,
     col, bg, font}{graphical parameters similar to standard
     \code{plot} functions; those applied only to points preceded with
     \code{p}, and those to text with \code{t}.}
  \item{\dots}{Other arguments passed to the functions.}

}

\details{

  \code{orditkplot} was moved to \CRANpkg{vegan3d} (version
  1.3-0). Install that package from CRAN and use in the old way.

}

\seealso{
  \code{\link{Deprecated}}
}
\keyword{misc}

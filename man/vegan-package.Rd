\name{vegan-package}
\alias{vegan-package}
\alias{vegan}
\docType{package}
\title{
The Community Ecology Package
}

\description{ The \pkg{vegan} package provides tools for descriptive
community ecology. It has most basic functions of diversity analysis,
community ordination and dissimilarity analysis. Most of its
multivariate tools can be used for other data types as well.  }

\details{The functions in the \pkg{vegan} package contain tools for
diversity analysis (see \code{\link{vignette}}
\code{vegandocs("diversity")}), ordination and analysis of
dissimilarities (see \code{\link{vignette}}
\code{vegandocs("intro")}). Together with the \pkg{labdsv} package,
the \pkg{vegan} package provides most standard tools of descriptive
community analysis. Package \pkg{ade4} provides an alternative
comprehensive package, and several other packages complement
\pkg{vegan} and provide tools for deeper analysis in specific
fields. Package \pkg{BiodiversityR} provides a GUI for a large subset
of \pkg{vegan} functionality.

The \pkg{vegan} package is developed at R-Forge
(\url{http://vegan.r-forge.r-project.org}).  The R-Forge provides
up-to-date information and mailing lists and fora for help queries and
bug reports. Bug reports can also be emailed to the function authors
of the package maintainers.

The \pkg{vegan} documents can be read with \code{\link{vegandocs}}
function. In addition to \code{\link{vignette}}s of basic usage, you
can also read \code{NEWS} on the new features and bug fixes in the
relese version (\code{vegandocs("NEWS")}), or more technical and fine
grained history in \code{vegandocs("ChangeLog")}.  The discussion on
design decisions can be read with \code{vegandocs("decision")}.  A
tutorial of the package at
\url{http://cc.oulu.fi/~jarioksa/opetus/metodi/vegantutor.pdf}
provides a more thorough introduction to the packaage.

To see the preferable citation of the package, type
\code{citation("vegan")}.  
}

\author{ The \pkg{vegan} development team is Jari Oksanen,
F. Guillaume Blanchet, Roeland Kindt, Pierre Legendre, R. B. O'Hara,
Gavin L. Simpson, Peter Solymos, M. Henry H. Stevens, Helene Wagner.

The maintainers at the R-Forge are Jari Oksanen <jari.oksanen@oulu.fi>
and Gavin Simpson <gavin.simpson@ucl.ac.uk>.
}

\keyword{ package }
\keyword{ multivariate }
\keyword{ univar }
\keyword{ models }
\keyword{ spatial }
\keyword{ nonparametric }
\keyword{ htest }
\keyword{ regression }
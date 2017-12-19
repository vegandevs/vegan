\name{model.matrix.cca}
\Rdversion{1.1}
\alias{model.matrix.cca}
\alias{model.matrix.rda}
\alias{model.frame.cca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Reconstruct Model Frame and Model Matrices of Constrained Ordination
}
\description{
  Function \code{model.frame.cca} reconstructs a \code{\link{data.frame}}
  with the variables used in the constrained ordination method 
  (\code{\link{cca}}, \code{\link{rda}} or \code{\link{capscale}}. 
  Function \code{model.matrix.cca} creates a list of design matrices used 
  in constrained ordination. The items of the list are called 
  \code{Conditions} and \code{Constraints}. If either partial 
  (\code{Conditions}) or constrained component was missing, 
  a single matrix is returned. 
}
\usage{
\method{model.frame}{cca}(formula, ...)
\method{model.matrix}{cca}(object, ...)
}

\arguments{
  \item{formula, object}{
  A constrained ordination result object from which the needed
  information is extracted.
}
  \item{\dots}{
  Other arguments passed to the default method of the function.
}
}
\details{
  The constrained ordination method objects do not save data on model 
  frame or model matrix, and the functions must reconstruct the
  information from the result object. This is usually possible for the
  model matrix. The reconstruction of model frame is only possible if
  the model was fitted with \code{formula} and \code{data} arguments,
  and still fails if the \code{data} are unavailable.
}
\value{
  Returns a data frame (\code{model.frame}) or an unnamed matrix or a 
  list of two matrices called \code{Constraints} and \code{Conditions} 
  (\code{model.matrix}). 
}

\note{The \code{model.matrix} returns the unweighted model matrix also
  for \code{\link{cca}}. Prior to \pkg{vegan} version 2.5-0 it returned
  the weighted model matrix.}

\author{
 Jari Oksanen
}
\seealso{
 \code{\link{model.frame}}, \code{\link{model.matrix}}. The model matrix
 is built from the QR decomposition component of \code{\link{cca.object}}
 using \code{\link{qr.X}}.
}
\examples{
data(dune)
data(dune.env)
mod <- cca(dune ~  poly(A1, 2) + Management + Use, dune.env)
model.frame(mod)
model.matrix(mod)
}
\keyword{ models }
\keyword{ multivariate }

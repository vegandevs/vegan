\name{MDSaddpoints}
\alias{MDSaddpoints}
\alias{dist2xy}
\title{Add New Points to NMDS ordination}

\description{Add new points to an existing \code{\link{metaMDS}} or
  \code{\link{monoMDS}} ordination.}

\usage{
MDSaddpoints(nmds, dis, neighbours = 5, maxit = 200)

dist2xy(dist, pick, type = c("xy", "xx"), invert = FALSE)
}
\arguments{

\item{nmds}{Result object from \code{\link{metaMDS}} or
  \code{\link{monoMDS}}. The configuration of points is fixed,
  but new points are added.}

\item{dis}{Rectangular non-symmetric dissimilarity matrix among new
  points (rows) and old fixed points (columns). Such matrix can be
  extracted from complete dissimilarities of both old and new points
  with \code{dist2xy}, or calculated with \code{\link{designdist2}}.}

\item{neighbours}{Number of nearest points used to get the starting
  locations for new points.}

\item{maxit}{Maximum number of iterations.}

\item{dist}{Input dissimilarities.}

\item{pick}{Indices (integers) of selected observations or a logical
  vector that is \code{TRUE} for picked items. The output will be
  in the original order and will not be reordered by this argument.}

\item{type}{\code{"xy"} returns rectangular data of picked against
  not picked observations, and \code{"xx"} a subset of symmetric
  dissimilarities.}

\item{invert}{Invert \code{pick}: drop elements listed.}
}
\value{
Function return a list of class \code{"nmds"} (there are no other
objects of that type in \pkg{vegan}) with following elements

  \item{points}{Coordinates of added new points}
  \item{seeds}{Starting coordinates for new points.}
  \item{deltastress}{Change of stress with added points.}
  \item{iters}{Number of iterations.}
  \item{cause}{Cause of termination of iterations. Integer for
    convergence criteria in \code{\link{monoMDS}}.}

}

\details{

Function provides an interface to \code{monoMDS} Fortran code to add new
points to an existing ordination that will be regarded as fixed. The
function has a similar role as \code{predict} functions with
\code{newdata} in Euclidean ordination (e.g. \code{\link{predict.cca}}).
Input data must be a rectangular matrix of distances among new added
points (rows) and all fixed old points (columns). Such matrices can be
extracted from complete dissimilarities with helper function
\code{dist2xy}. Function \code{\link{designdist2}} can directly
calculate such rectangular dissimilarity matrices between sampling units
(rows) in two matries. In addition, \CRANpkg{analogue} has distance
function that can calculate dissimilarities among two matrices,
including functions that cannot be specified in
\code{\link{designdist2}}.

Great care is needed in preparing the dissimilarities for the
input. The dissimilarity index must be exactly the same as in the
fixed ordination, and columns must match old fixed points, and rows
added new points.

}
\examples{
## Cross-validation: remove a point when performing NMDS and add as
## a new points
data(dune)
d <- vegdist(dune)
## remove point 3 from ordination
mod3 <- metaMDS(dist2xy(d, 3, "xx", invert = TRUE), trace=0)
## add point 3 to the result
MDSaddpoints(mod3, dist2xy(d, 3))
## Use designdist2
d15 <- designdist(dune[1:15,])
m15 <- metaMDS(d15, trace=0)
MDSaddpoints(m15, designdist2(dune[1:15,], dune[16:20,]))
}

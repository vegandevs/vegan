\encoding{UTF-8}
\name{avgdist}
\alias{avgdist}
\title{Averaged Subsampled Dissimilarity Matrices}
\description{
  The function computes the dissimilarity matrix of a dataset multiple
  times using \code{\link{vegdist}} while randomly subsampling the
  dataset each time. All of the subsampled iterations are then averaged
  (mean) to provide a distance matrix that represents the average of
  multiple subsampling iterations. This emulates the behavior of the
  distance matrix calculator within the Mothur microbial ecology toolkit.
}

\usage{avgdist(x, depth, distfun = vegdist, iterations = 100, dmethod = "bray", ...) } 
\arguments{
  \item{x}{Community data matrix.}
  \item{depth}{The subsampling depth to be used in each iteration.}
  \item{distfun}{The dissimilarity matrix function to be used. Default is the vegan \code{\link{vegdist}}}
  \item{iterations}{The number of random iterations to perform before averaging. Default is 100 iterations.}
  \item{dmethod}{Dissimilarity index to be used with the specified dissimilarity matrix function.}
}

\author{ Geoffrey Hannigan. }

\note{
  The function builds on the function \code{\link{rrarefy}} and and additional
  distance matrix function (e.g. \code{\link{vegdist}}) to add more meaningful representations of
  distances among randomly subsampled datasets by presenting the
  average of multiple random iterations. This function runs using
  the \code{\link{vegdist}} options of `diag = TRUE` and `upper =
  TRUE`. This functionality has been utilized in the Mothur standalone
  microbial ecology toolkit \href{https://www.mothur.org/wiki/Dist.shared#subsample}{here}.
}

\seealso{
 This function utilizes the \code{\link{vegdist}} and \code{\link{rrarefy}}
 vegan functions.}

\examples{
data(varespec)
# Multiply decimals to remove decimals and simulate count data for this example
countdist <- varespec * 1000
mean.avg.dist <- avgdist(countdist, depth = 500, iterations = 10)
mean.avg.dist
}
\keyword{ distance }

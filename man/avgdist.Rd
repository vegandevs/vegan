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

\usage{avgdist(x, depth, iterations = 100, dmethod = "bray") } 
\arguments{
  \item{x}{Community data matrix.}
  \item{dmethod}{Dissimilarity index, partial match to
    \code{"manhattan"}, \code{"euclidean"}, \code{"canberra"},
    \code{"clark"}, \code{"bray"}, \code{"kulczynski"},
    \code{"jaccard"}, \code{"gower"}, \code{"altGower"},
    \code{"morisita"}, \code{"horn"}, \code{"mountford"}, \code{"raup"},
    \code{"binomial"}, \code{"chao"}, \code{"cao"} or \code{"mahalanobis"}.}
  \item{depth}{The subsampling depth to be used in each iteration.}
  \item{iterations}{The number of random iterations to perform before averaging.}
}

\author{ Geoffrey Hannigan. }

\note{
  The function builds on the functions \code{\link{vegdist}} and
  \code{\link{rrarefy}} to add more meaningful representations of
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

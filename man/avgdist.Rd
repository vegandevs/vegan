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

\usage{avgdist(x, sample, distfun = vegdist, meanfun = mean, transf = NULL, iterations = 100, dmethod = "bray", ...) } 
\arguments{
  \item{x}{Community data matrix.}
  \item{sample}{The subsampling depth to be used in each iteration.}
  \item{distfun}{The dissimilarity matrix function to be used. Default is the vegan \code{\link{vegdist}}}
  \item{meanfun}{The calculation to use for the average (mean or median).}
  \item{transf}{Option for transforming the count data before calculating the distance matrix. Any base transformation option can be used (e.g. \code{\link{sqrt}})}
  \item{iterations}{The number of random iterations to perform before averaging. Default is 100 iterations.}
  \item{dmethod}{Dissimilarity index to be used with the specified dissimilarity matrix function.}
  \item{...}{Any additional arguments to add to the distance function specified.}
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
# Import an example count dataset
data(BCI)
# Test the base functionality
mean.avg.dist <- avgdist(BCI, sample = 50, iterations = 10)
# Test the transformation function
mean.avg.dist.t <- avgdist(BCI, sample = 50, iterations = 10, transf = sqrt)
# Test the median functionality
median.avg.dist <- avgdist(BCI, sample = 50, iterations = 10, meanfun = median)
# Print the resulting tables
mean.avg.dist
mean.avg.dist.t
median.avg.dist
}
\keyword{ distance }

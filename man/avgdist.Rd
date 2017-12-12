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

\usage{
avgdist(x, sample, distfun = vegdist, meanfun = mean,
    transf = NULL, iterations = 100, dmethod = "bray", ...)
}
\arguments{
  \item{x}{Community data matrix.}
  \item{sample}{The subsampling depth to be used in each iteration. Samples that do not meet this threshold will be removed from the analysis, and their identify returned to the user in stdout.}
  \item{distfun}{The dissimilarity matrix function to be used. Default is the vegan \code{\link{vegdist}}}
  \item{meanfun}{The calculation to use for the average (mean or median).}
  \item{transf}{Option for transforming the count data before calculating the distance matrix. Any base transformation option can be used (e.g. \code{\link{sqrt}})}
  \item{iterations}{The number of random iterations to perform before averaging. Default is 100 iterations.}
  \item{dmethod}{Dissimilarity index to be used with the specified dissimilarity matrix function. Default is Bray-Curtis}
  \item{...}{Any additional arguments to add to the distance function or mean/median function specified.}
}

\author{ Geoffrey Hannigan. }

\note{

  The function builds on the function \code{\link{rrarefy}} and and
  additional distance matrix function (e.g. \code{\link{vegdist}}) to
  add more meaningful representations of distances among randomly
  subsampled datasets by presenting the average of multiple random
  iterations. This function runs using the \code{\link{vegdist}}. This
  functionality has been utilized in the Mothur standalone microbial
  ecology toolkit
  \href{https://www.mothur.org/wiki/Dist.shared#subsample}{here}.

}

\seealso{
 This function utilizes the \code{\link{vegdist}} and \code{\link{rrarefy}}
 functions.}

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
head(as.matrix(mean.avg.dist))
head(as.matrix(mean.avg.dist.t))
head(as.matrix(median.avg.dist))
# Run example to illustrate low variance of mean, median, and stdev results
# Mean and median std dev are around 0.05
sdd <- avgdist(BCI, sample = 50, iterations = 100, meanfun = sd)
summary(mean.avg.dist)
summary(median.avg.dist)
summary(sdd)
# Test for when subsampling depth excludes some samples
# Return samples that are removed for not meeting depth filter
depth.avg.dist <- avgdist(BCI, sample = 450, iterations = 10)
# Print the result
depth.avg.dist
}
\keyword{ multivariate }

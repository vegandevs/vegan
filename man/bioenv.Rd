\name{bioenv}
\alias{bioenv}
\alias{bioenv.default}
\alias{bioenv.formula}
\alias{summary.bioenv}
\alias{bioenvdist}

\title{Best Subset of Environmental Variables with
  Maximum (Rank) Correlation with Community Dissimilarities }
\description{
  Function finds the best subset of environmental variables, so that
  the Euclidean distances of scaled environmental variables have the
  maximum (rank) correlation with community dissimilarities.  
}
\usage{
\method{bioenv}{default}(comm, env, method = "spearman", index = "bray",
       upto = ncol(env), trace = FALSE, partial = NULL, 
       metric = c("euclidean", "mahalanobis", "manhattan", "gower"),
       parallel = getOption("mc.cores"), ...)
\method{bioenv}{formula}(formula, data, ...)
bioenvdist(x, which = "best")
}

\arguments{
  \item{comm}{Community data frame or a dissimilarity object or a square
    matrix that can be interpreted as dissimilarities. }
  \item{env}{Data frame of continuous environmental variables. }
  \item{method}{The correlation method used in \code{\link{cor}}.}
  \item{index}{The dissimilarity index used for community data (\code{comm}) 
    in \code{\link{vegdist}}. This is ignored if \code{comm} are dissimilarities.}
  \item{upto}{Maximum number of parameters in studied subsets.}
  \item{formula, data}{Model \code{\link{formula}} and data.}
  \item{trace}{Trace the calculations }
  \item{partial}{Dissimilarities partialled out when inspecting
    variables in \code{env}.}
  \item{metric}{Metric used for distances of environmental distances. See 
    Details.}
  \item{parallel}{Number of parallel processes or a predefined socket
    cluster.  With \code{parallel = 1} uses ordinary, non-parallel
    processing. The parallel processing is done with \pkg{parallel}
    package.}
  \item{x}{\code{bioenv} result object.}
  \item{which}{The number of the model for which the environmental
    distances are evaluated, or the \code{"best"} model.}
  \item{...}{Other arguments passed to \code{\link{cor}}.}
}
\details{
  
  The function calculates a community dissimilarity matrix using
  \code{\link{vegdist}}.  Then it selects all possible subsets of
  environmental variables, \code{\link{scale}}s the variables, and
  calculates Euclidean distances for this subset using
  \code{\link{dist}}.  The function finds the correlation between
  community dissimilarities and environmental distances, and for each
  size of subsets, saves the best result.  There are \eqn{2^p-1}
  subsets of \eqn{p} variables, and an exhaustive search may take a
  very, very, very long time (parameter \code{upto} offers a partial
  relief).

  The argument \code{metric} defines distances in the given set of
  environmental variables.  With \code{metric = "euclidean"}, the
  variables are scaled to unit variance and Euclidean distances are
  calculated. With \code{metric = "mahalanobis"}, the Mahalanobis
  distances are calculated: in addition to scaling to unit variance,
  the matrix of the current set of environmental variables is also
  made orthogonal (uncorrelated). With \code{metric = "manhanttan"},
  the variables are scaled to unit range and Manhattan distances are
  calculated, so that the distances are sums of differences of
  environmental variables.  With \code{metric = "gower"}, the Gower
  distances are calculated using function
  \code{\link[cluster]{daisy}}. This allows also using factor
  variables, but with continuous variables the results are equal to
  \code{metric = "manhattan"}.

  The function can be called with a model \code{\link{formula}} where
  the LHS is the data matrix and RHS lists the environmental variables.
  The formula interface is practical in selecting or transforming
  environmental variables.

  With argument \code{partial} you can perform \dQuote{partial}
  analysis. The partializing item must be a dissimilarity object of
  class \code{\link{dist}}. The
  \code{partial} item can be used with any correlation \code{method},
  but it is strictly correct only for Pearson.

  Function \code{bioenvdist} recalculates the environmental distances
  used within the function. The default is to calculate distances for
  the best model, but the number of any model can be given.
  
  Clarke & Ainsworth (1993) suggested this method to be used for
  selecting the best subset of environmental variables in interpreting
  results of nonmetric multidimensional scaling (NMDS). They recommended a
  parallel display of NMDS of community dissimilarities and NMDS of
  Euclidean distances from the best subset of scaled environmental
  variables.  They warned against the use of Procrustes analysis, but
  to me this looks like a good way of comparing these two ordinations.

  Clarke & Ainsworth wrote a computer program BIO-ENV giving the name to
  the current function. Presumably BIO-ENV
  was later incorporated in Clarke's PRIMER software (available for
  Windows).  In addition, Clarke & Ainsworth suggested a novel method of
  rank correlation which is not available in the current function.
}

\value{
  The function returns an object of class \code{bioenv} with a
  \code{summary} method.
}

\references{
  Clarke, K. R & Ainsworth, M. 1993. A method of linking multivariate
  community structure to environmental variables. \emph{Marine Ecology
    Progress Series}, 92, 205--219.
}
\author{ Jari Oksanen }

\note{ If you want to study the \sQuote{significance} of \code{bioenv}
  results, you can use function \code{\link{mantel}} or
  \code{\link{mantel.partial}} which use the same definition of
  correlation.  However, \code{bioenv} standardizes environmental
  variables depending on the used metric, and you must do the same in
  \code{\link{mantel}} for comparable results (the standardized data are
  returned as item \code{x} in the result object). It is safest to use
  \code{bioenvdist} to extract the environmental distances that really
  were used within \code{bioenv}. NB., \code{bioenv} selects variables
  to maximize the Mantel correlation, and significance tests based on
  \emph{a priori} selection of variables are biased.  }

\seealso{\code{\link{vegdist}}, \code{\link{dist}}, \code{\link{cor}}
  for underlying routines, \code{\link{monoMDS}} and
  \code{\link{metaMDS}} for ordination, \code{\link{procrustes}} for
  Procrustes analysis, \code{\link{protest}} for an alternative, and
  \code{\link{rankindex}} for studying alternatives to the default
  Bray-Curtis index.}

\examples{
# The method is very slow for large number of possible subsets.
# Therefore only 6 variables in this example.
data(varespec)
data(varechem)
sol <- bioenv(wisconsin(varespec) ~ log(N) + P + K + Ca + pH + Al, varechem)
sol
summary(sol)
}
\keyword{ multivariate }


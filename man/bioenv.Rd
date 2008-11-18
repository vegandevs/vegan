\name{bioenv}
\alias{bioenv}
\alias{bioenv.default}
\alias{bioenv.formula}
\alias{print.bioenv}
\alias{summary.bioenv}
\alias{print.summary.bioenv}
\alias{ripley.subsets}
\alias{ripley.subs}

\title{Best Subset of Environmental Variables with
  Maximum (Rank) Correlation with Community Dissimilarities }
\description{
  Function finds the best subset of environmental variables, so that
  the Euclidean distances of scaled environmental variables have the
  maximum (rank) correlation with community dissimilarities.  
}
\usage{
\method{bioenv}{default}(comm, env, method = "spearman", index = "bray",
       upto = ncol(env), trace = FALSE, partial = NULL, ...)
\method{bioenv}{formula}(formula, data, ...)
}

\arguments{
  \item{comm}{Community data frame. }
  \item{env}{Data frame of continuous environmental variables. }
  \item{method}{The correlation method used in \code{\link{cor}}.}
  \item{index}{The dissimilarity index used for community data in
    \code{\link{vegdist}}. }
  \item{upto}{Maximum number of parameters in studied subsets.}
  \item{formula, data}{Model \code{\link{formula}} and data.}
  \item{trace}{Trace the advance of calculations }
  \item{partial}{Dissimilarities partialled out when inspecting
    variables in \code{env}.}
  \item{...}{Other arguments passed to \code{\link{cor}}.}
}
\details{
  The function calculates a community dissimilarity matrix using
  \code{\link{vegdist}}.  Then it selects all possible subsets of
  environmental variables, \code{\link{scale}}s the variables, and
  calculates Euclidean distances for this subset using
  \code{\link{dist}}.  Then it finds the correlation between community
  dissimilarities and environmental distances, and for each size of
  subsets, saves the best result. 
  There are \eqn{2^p-1} subsets of \eqn{p} variables, and an exhaustive
  search may take a very, very, very long time (parameter \code{upto} offers a
  partial relief). 

  The function can be called with a model \code{\link{formula}} where
  the LHS is the data matrix and RHS lists the environmental variables.
  The formula interface is practical in selecting or transforming
  environmental variables.

  With argument \code{partial} you can perform \dQuote{partial}
  analysis. The partializing item must be a dissimilarity object of
  class \code{\link{dist}}. The
  \code{partial} item can be used with any correlation \code{method},
  but it is strictly correct only for Pearson.
  
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
\author{ Jari Oksanen. The code for selecting all possible subsets was
  posted to the R mailing list by Prof. B. D. Ripley in 1999. }

\note{
  If you want to study the \sQuote{significance} of \code{bioenv}
  results, you can use function \code{\link{mantel}} or
  \code{\link{mantel.partial}} which use the same definition of
  correlation. 
  However, \code{bioenv} standardizes environmental variables to unit standard
  deviation using function \code{\link{scale}} and you must do the same
  in \code{\link{mantel}} for comparable results. Further, \code{bioenv}
  selects variables to maximize the Mantel correlation, and significance
  tests based on \emph{a priori} selection of variables are biased. 
  }

\seealso{\code{\link{vegdist}},
  \code{\link{dist}}, \code{\link{cor}} for underlying routines,
  \code{\link[MASS]{isoMDS}} for ordination, \code{\link{procrustes}}
  for Procrustes analysis, \code{\link{protest}} for an alternative, and
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


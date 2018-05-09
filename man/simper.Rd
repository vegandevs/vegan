\encoding{UTF-8}
\name{simper}
\alias{simper}
\alias{summary.simper}
\title{Similarity Percentages}

\description{
  Discriminating species between two groups using
  Bray-Curtis dissimilarities
}

\usage{
simper(comm, group, permutations = 0, trace = FALSE, 
    parallel = getOption("mc.cores"), ...)
\method{summary}{simper}(object, ordered = TRUE,
    digits = max(3,getOption("digits") - 3), ...)
}

\arguments{
  \item{comm}{Community data matrix.}
  \item{group}{Factor describing the group structure. Must have at
    least 2 levels.}
  \item{permutations}{a list of control values for the permutations
    as returned by the function \code{\link[permute]{how}}, or the
    number of permutations required, or a permutation matrix where each
    row gives the permuted indices.}
  \item{trace}{Trace permutations.}
  \item{object}{an object returned by \code{simper}.}
  \item{ordered}{Logical; Should the species be ordered by their
    average contribution?}
  \item{digits}{Number of digits in output.}
  \item{parallel}{Number of parallel processes or a predefined socket
    cluster.  With \code{parallel = 1} uses ordinary, non-parallel
    processing.}
  \item{...}{Parameters passed to other functions. In \code{simper} the
    extra parameters are passed to \code{\link[permute]{shuffleSet}} if
    permutations are used.}
}

\details{ Similarity percentage, \code{simper} (Clarke 1993) is based
  on the decomposition of Bray-Curtis dissimilarity index (see
  \code{\link{vegdist}}, \code{\link{designdist}}). The contribution
  of individual species \eqn{i} to the overall Bray-Curtis dissimilarity
  \eqn{d_{jk}}{d[jk]} is given by

  \deqn{d_{ijk} = \frac{|x_{ij}-x_{ik}|}{\sum_{i=1}^S (x_{ij}+x_{ik})}}{d[ijk] = abs(x[ij]-x[ik])/sum(x[ij]+x[ik])}
  
  where \eqn{x} is the abundance of species \eqn{i} in sampling units
  \eqn{j} and \eqn{k}. The overall index is the sum of the individual
  contributions over all \eqn{S} species 
  \eqn{d_{jk}=\sum_{i=1}^S d_{ijk}}{d[jk] = sum(i=1..S) d[ijk]}. 
  
  The \code{simper} functions performs pairwise comparisons of groups
  of sampling units and finds the average contributions
  of each species to the average overall Bray-Curtis dissimilarity.

  The function displays most important species for each pair of
  \code{groups}.  These species contribute at least to 70 \% of the
  differences between groups.  The function returns much more
  extensive results which can be accessed directly from the result
  object (see section Value). Function \code{summary} transforms the
  result to a list of data frames. With argument \code{ordered = TRUE}
  the data frames also include the cumulative contributions and
  are ordered by species contribution.

  The results of \code{simper} can be very difficult to interpret. The
  method very badly confounds the mean between group differences and
  within group variation, and seems to single out variable species
  instead of distinctive species (Warton et al. 2012). Even if you make
  groups that are copies of each other, the method will single out
  species with high contribution, but these are not contributions
  to non-existing between-group differences but to within-group
  variation in species abundance.
  

}

\value{
  A list of class \code{"simper"} with following items:
  \item{species}{The species names.}
  \item{average}{Average contribution to overall dissimilarity.}
  \item{overall}{The overall between-group dissimilarity.} 
  \item{sd}{Standard deviation of contribution.} 
  \item{ratio}{Average to sd ratio.}
  \item{ava, avb}{Average abundances per group.}
  \item{ord}{An index vector to order vectors by their contribution or
    order \code{cusum} back to the original data order.}
  \item{cusum}{Ordered cumulative contribution.}
  \item{p}{Permutation \eqn{p}-value. Probability of getting a larger
    or equal average contribution in random permutation of the group
    factor. These area only available if \code{permutations} were used
    (default: not calculated).}

}

\examples{
data(dune)
data(dune.env)
(sim <- with(dune.env, simper(dune, Management)))
summary(sim)
}
\author{
  Eduard Szöcs \email{eduardszoecs@gmail.com}
}

\references{
  Clarke, K.R. 1993. Non-parametric multivariate analyses of changes
    in community structure. \emph{Australian Journal of Ecology}, 18,
    117–143.

  Warton, D.I., Wright, T.W., Wang, Y. 2012. Distance-based multivariate
    analyses confound location and dispersion effects. \emph{Methods in
    Ecology and Evolution}, 3, 89--101.
}
\keyword{multivariate}


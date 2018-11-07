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
simper(comm, group, permutations = 999, parallel = 1, ...)
\method{summary}{simper}(object, ordered = TRUE,
    digits = max(3,getOption("digits") - 3), ...)
}

\arguments{
  \item{comm}{Community data.}
  \item{group}{Factor describing the group structure. If this is
    missing or has only one level, contributions are estimated for
    non-grouped data and dissimilarities only show the overall
    heterogeneity in species abundances.}
  \item{permutations}{a list of control values for the permutations
    as returned by the function \code{\link[permute]{how}}, or the
    number of permutations required, or a permutation matrix where each
    row gives the permuted indices.}
  \item{object}{an object returned by \code{simper}.}
  \item{ordered}{Logical; Should the species be ordered by their
    average contribution?}
  \item{digits}{Number of digits in output.}
  \item{parallel}{Number of parallel processes or a predefined socket
    cluster.  With \code{parallel = 1} uses ordinary, non-parallel
    processing. (Not yet implemented).}
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
  of sampling units and finds the contribution of each species to the
  average between-group Bray-Curtis dissimilarity. Although the method
  is called \dQuote{Similarity Percentages}, it really studied
  dissimilarities instead of similarities (Clarke 1993).

  The function displays most important species for each pair of
  \code{groups}.  These species contribute at least to 70 \% of the
  differences between groups.  The function returns much more
  extensive results (including all species) which can be accessed
  directly from the result object (see section Value). Function
  \code{summary} transforms the result to a list of data frames. With
  argument \code{ordered = TRUE} the data frames also include the
  cumulative contributions and are ordered by species contribution.

  The results of \code{simper} can be very difficult to interpret and
  they are often misunderstood even in publications. The method gives
  the contribution of each species to overall dissimilarities, but
  these are caused by variation in species abundances, and only partly
  by differences among groups.  Even if you make groups that are
  copies of each other, the method will single out species with high
  contribution, but these are not contributions to non-existing
  between-group differences but to random noise variation in species
  abundances. The most abundant species usually have highest
  variances, and they have high contributions even when they do not
  differ among groups. Permutation tests study the differences among
  groups, and they can be used to find out the species for which the
  differences among groups is an important component of their
  contribution to dissimilarities. Analysis without \code{group}
  argument will find species contributions to the average overall
  dissimilarity among sampling units. These non-grouped contributions
  can be compared to grouped contributions to see how much added value
  the grouping has for each species.

}

\value{
  A list of class \code{"simper"} with following items:
  \item{species}{The species names.}
  \item{average}{Species contribution to average between-group dissimilarity.}
  \item{overall}{The average between-group dissimilarity. This is the sum of
    the item \code{average}.}
  \item{sd}{Standard deviation of contribution.} 
  \item{ratio}{Average to sd ratio.}
  \item{ava, avb}{Average abundances per group.}
  \item{ord}{An index vector to order vectors by their contribution or
    order \code{cusum} back to the original data order.}
  \item{cusum}{Ordered cumulative contribution. These are based on item
    \code{average}, but they sum up to total 1.}
  \item{p}{Permutation \eqn{p}-value. Probability of getting a larger
    or equal average contribution in random permutation of the group
    factor. These area only available if \code{permutations} were used
    (default: not calculated).}
}

\seealso{Function \code{\link{meandist}} shows the average between-group
   dissimilarities (as well as the within-group dissimilarities).}

\examples{
data(dune)
data(dune.env)
(sim <- with(dune.env, simper(dune, Management, permutations = 99)))
summary(sim)
}
\author{
  Eduard Szöcs and Jari Oksanen.
}

\references{
  Clarke, K.R. 1993. Non-parametric multivariate analyses of changes
    in community structure. \emph{Australian Journal of Ecology}, 18,
    117–143.
}
\keyword{multivariate}


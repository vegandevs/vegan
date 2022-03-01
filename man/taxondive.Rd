\name{taxondive}
\alias{taxondive}
\alias{summary.taxondive}
\alias{plot.taxondive}
\alias{taxa2dist}

\title{ Indices of Taxonomic Diversity and Distinctness }
\description{
  Function finds indices of taxonomic diversity and distinctness,
  which are averaged taxonomic distances among species or individuals in
  the community (Clarke & Warwick 1998, 2001)
}
\usage{
taxondive(comm, dis, match.force = FALSE)
taxa2dist(x, varstep = FALSE, check = TRUE, labels)
}

\arguments{
  \item{comm}{Community data. }
  \item{dis}{Taxonomic distances among taxa in \code{comm}. This should
    be a \code{\link{dist}} object or a symmetric square matrix. }
  \item{match.force}{Force matching of column names in \code{comm} and
    labels in \code{dis}. If \code{FALSE}, matching only happens when
    dimensions differ, and in that case the species must be in identical
    order in both.}
  \item{x}{Classification table with a row for each species or other
    basic taxon, and columns
    for identifiers of its classification at higher levels.}
  \item{varstep}{Vary step lengths between successive levels
    relative to proportional loss of the number of distinct classes.}
  \item{check}{If \code{TRUE}, remove all redundant levels which are
    different for all rows or constant for all rows and regard each row
    as a different basal taxon (species). If \code{FALSE} all
    levels are retained and basal taxa (species) also must be coded as
    variables (columns). You will get a warning if species are not
    coded, but you can ignore this if that was your intention.}
  \item{labels}{The \code{labels} attribute of taxonomic distances. Row
    names will be used if this is not given. Species will be matched by
    these labels in \code{comm} and \code{dis} in \code{taxondive} if
    these have different dimensions.}
}
\details{
  Clarke & Warwick (1998, 2001) suggested several alternative indices of
  taxonomic diversity or distinctness. Two basic indices are called
  taxonomic diversity (\eqn{\Delta}) and distinctness (\eqn{\Delta^*}):
  \tabular{c}{
    \eqn{\Delta = (\sum \sum_{i<j} \omega_{ij} x_i x_j)/(n (n-1) / 2)} \cr
    \eqn{\Delta^* = (\sum \sum_{i<j} \omega_{ij} x_i x_j)/(\sum \sum_{i<j} x_i x_j) }
  }
  The equations give the index value for a single site, and summation
  goes over species \eqn{i} and \eqn{j}. Here \eqn{\omega} are taxonomic
  distances among taxa, and \eqn{x} are species abundances, and \eqn{n}
  is the total abundance for a site. 
  With presence/absence  data both indices reduce to the same index
  \eqn{\Delta^+}, and for this index Clarke & Warwick (1998) also have
  an estimate of its standard deviation. Clarke & Warwick (2001) 
  presented two new indices: \eqn{s\Delta^+} is the product of species
  richness and \eqn{\Delta^+}, and index of variation in
  taxonomic distinctness (\eqn{\Lambda^+}) defined as
  \tabular{c}{
    \eqn{\Lambda^+ = (\sum \sum_{i<j} \omega_{ij}^2)/(n (n-1) / 2) - (\Delta^+)^2}
  }

  The \code{dis} argument must be species dissimilarities. These must be
  similar to dissimilarities produced by \code{\link{dist}}. It is
  customary to have integer steps of taxonomic hierarchies, but other
  kind of dissimilarities can be used, such as those from phylogenetic
  trees or genetic differences.  Further, the \code{dis} need not be
  taxonomic, but other species classifications can be used. 

  Function \code{taxa2dist} can produce a suitable \code{dist} object
  from a classification table. Each species (or basic taxon) corresponds
  to a row of the classification table, and columns give the
  classification at different levels. With \code{varstep = FALSE} the
  successive levels will be separated by equal steps, and with
  \code{varstep = TRUE} the step length is relative to the proportional
  decrease in the number of classes (Clarke & Warwick 1999).
  With \code{check = TRUE}, the function removes classes which are distinct for all
  species or which combine all species into one class, and assumes that
  each row presents a distinct basic taxon. The function scales
  the distances so that longest path length between
  taxa is 100 (not necessarily when \code{check = FALSE}). 
  
  Function \code{plot.taxondive} plots \eqn{\Delta^+} against Number of
  species, together with expectation and its approximate 2*sd
  limits. Function \code{summary.taxondive} finds the \eqn{z} values and
  their significances from Normal distribution for \eqn{\Delta^+}.
}
\value{
Function returns an object of class \code{taxondive} with following items:
  \item{Species }{Number of species for each site.}
  \item{D, Dstar, Dplus, SDplus, Lambda}{\eqn{\Delta}, \eqn{\Delta^*},
    \eqn{\Delta^+},  \eqn{s\Delta^+} and \eqn{\Lambda^+} 
    for each site.}
  \item{sd.Dplus}{Standard deviation of \eqn{\Delta^+}.}
  \item{ED, EDstar, EDplus}{Expected values of corresponding
    statistics.}

  Function \code{taxa2dist} returns an object of class \code{"dist"}, with
  an attribute \code{"steps"} for the step lengths between successive levels.
}
\references{
 Clarke, K.R & Warwick, R.M. (1998) A taxonomic distinctness index and
 its statistical properties. \emph{Journal of Applied Ecology} 35,
 523--531.

 Clarke, K.R. & Warwick, R.M. (1999) The taxonomic distinctness measure
 of biodiversity: weighting of step lengths between hierarchical
 levels. \emph{Marine Ecology Progress Series} 184: 21--29.

 Clarke, K.R. & Warwick, R.M. (2001) A further biodiversity index
 applicable to species lists: variation in taxonomic
 distinctness. \emph{Marine Ecology Progress Series} 216, 265--278.
}
\author{ Jari Oksanen }
\note{
  The function is still preliminary and may change. The scaling of
  taxonomic dissimilarities influences the results. If you multiply
  taxonomic distances (or step lengths) by a constant, the values of all
  Deltas will be multiplied with the same constant, and the value of
  \eqn{\Lambda^+} by the square of the constant.
}
\seealso{\code{\link{diversity}}. }
\examples{
## Preliminary: needs better data and some support functions
data(dune)
data(dune.taxon)
# Taxonomic distances from a classification table with variable step lengths.
taxdis <- taxa2dist(dune.taxon, varstep=TRUE)
plot(hclust(taxdis), hang = -1)
# Indices
mod <- taxondive(dune, taxdis)
mod
summary(mod)
plot(mod)
}

\keyword{ univar }


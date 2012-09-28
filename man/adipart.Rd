\encoding{UTF-8}
\name{adipart}
\alias{adipart}
\alias{adipart.default}
\alias{adipart.formula}
\alias{hiersimu}
\alias{hiersimu.default}
\alias{hiersimu.formula}

\title{Additive Diversity Partitioning and Hierarchical Null Model Testing}
\description{
In additive diversity partitioning, mean values of alpha diversity at lower levels of a sampling 
hierarchy are compared to the total diversity in the entire data set (gamma diversity). 
In hierarchical null model testing, a statistic returned by a function is evaluated 
according to a nested hierarchical sampling design (\code{hiersimu}).
}
\usage{
adipart(...)
\method{adipart}{default}(y, x, index=c("richness", "shannon", "simpson"),
    weights=c("unif", "prop"), relative = FALSE, nsimul=99, ...)
\method{adipart}{formula}(formula, data, index=c("richness", "shannon", "simpson"),
    weights=c("unif", "prop"), relative = FALSE, nsimul=99, ...)

hiersimu(...)
\method{hiersimu}{default}(y, x, FUN, location = c("mean", "median"),
    relative = FALSE, drop.highest = FALSE, nsimul=99, ...)
\method{hiersimu}{formula}(formula, data, FUN, location = c("mean", "median"),
    relative = FALSE, drop.highest = FALSE, nsimul=99, ...)
}
\arguments{
  \item{y}{A community matrix.}
  \item{x}{A matrix with same number of rows as in \code{y}, columns
    coding the levels of sampling hierarchy. The number of groups within
    the hierarchy must decrease from left to right. If \code{x} is missing,
    function performs an overall decomposition into alpha, beta and
    gamma diversities.}
  \item{formula}{A two sided model formula in the form \code{y ~ x},
    where \code{y} is the community data matrix with samples as rows and
    species as column. Right hand side (\code{x}) must grouping vaiables
    referring to levels of sampling hierarchy, terms from right to left
    will be treated as nested (first column is the lowest, last is the
    highest level, at least two levels specified). Interaction terms are
    not allowed.}

  \item{data}{A data frame where to look for variables defined in the
    right hand side of \code{formula}. If missing, variables are looked
    in the global environment.}

  \item{index}{Character, the diversity index to be calculated (see Details).}

  \item{weights}{Character, \code{"unif"} for uniform weights,
    \code{"prop"} for weighting proportional to sample abundances to use
    in weighted averaging of individual alpha values within strata of a
    given level of the sampling hierarchy.}

  \item{relative}{Logical, if \code{TRUE} then alpha and beta diversity
    values are given relative to the value of gamma for function
    \code{adipart}.}

  \item{nsimul}{Number of permutations to use if \code{matr} is not of
    class 'permat'.  If \code{nsimul = 0}, only the \code{FUN} argument
    is evaluated. It is thus possible to reuse the statistic values
    without using a null model.}

  \item{FUN}{A function to be used by \code{hiersimu}. This must be
    fully specified, because currently other arguments cannot be passed
    to this function via \code{\dots}.}

  \item{location}{Character, identifies which function (mean or median)
    is to be used to calculate location of the samples.}

  \item{drop.highest}{Logical, to drop the highest level or not. When
    \code{FUN} evaluates only arrays with at least 2 dimensions, highest
    level should be dropped, or not selected at all.}

  \item{\dots}{Other arguments passed to functions, e.g. base of
    logarithm for Shannon diversity, or \code{method}, \code{thin} or
    \code{burnin} arguments for \code{\link{oecosimu}}.}
}

\details{

  Additive diversity partitioning means that mean alpha and beta
  diversities add up to gamma diversity, thus beta diversity is measured
  in the same dimensions as alpha and gamma (Lande 1996). This additive
  procedure is then extended across multiple scales in a hierarchical
  sampling design with \eqn{i = 1, 2, 3, \ldots, m} levels of sampling
  (Crist et al. 2003). Samples in lower hierarchical levels are nested
  within higher level units, thus from \eqn{i=1} to \eqn{i=m} grain size
  is increasing under constant survey extent. At each level \eqn{i},
  \eqn{\alpha_i} denotes average diversity found within samples.

  At the highest sampling level, the diversity components are calculated
  as \deqn{\beta_m  = \gamma -  \alpha_m}{beta_m = gamma -  alpha_m} For
  each  lower   sampling  level   as  \deqn{\beta_i  =   \alpha_{i+1}  -
  \alpha_i}{beta_i =  alpha_i+1 - alpha_i} Then,  the additive partition
  of diversity is \deqn{\gamma  = \alpha_1 + \sum_{i=1}^m \beta_i}{gamma
  = alpha_1 + sum(beta_i)}

  Average alpha components can be weighted uniformly
  (\code{weight="unif"}) to calculate it as simple average, or
  proportionally to sample abundances (\code{weight="prop"}) to
  calculate it as weighted average as follows \deqn{\alpha_i =
  \sum_{j=1}^{n_i} D_{ij} w_{ij}}{alpha_i = sum(D_ij*w_ij)} where
  \eqn{D_{ij}} is the diversity index and \eqn{w_{ij}} is the weight
  calculated for the \eqn{j}th sample at the \eqn{i}th sampling level.

  The implementation of additive diversity partitioning in
  \code{adipart} follows Crist et al. 2003. It is based on species
  richness (\eqn{S}, not \eqn{S-1}), Shannon's and Simpson's diversity
  indices stated as the \code{index} argument.

  The expected diversity components are calculated \code{nsimul} times
  by individual based randomisation of the community data matrix. This
  is done by the \code{"r2dtable"} method in \code{\link{oecosimu}} by
  default.

  \code{hiersimu} works almost in the same way as \code{adipart}, but
  without comparing the actual statistic values returned by \code{FUN}
  to the highest possible value (cf. gamma diversity).  This is so,
  because in most of the cases, it is difficult to ensure additive
  properties of the mean statistic values along the hierarchy.

}
\value{

  An object of class \code{"adipart"} or \code{"hiersimu"} with same
  structure as \code{\link{oecosimu}} objects.

}

\references{

  Crist,   T.O.,   Veech,    J.A.,   Gering,   J.C.   and   Summerville,
  K.S.  (2003).  Partitioning species  diversity  across landscapes  and
  regions:  a hierarchical  analysis of  \eqn{\alpha},  \eqn{\beta}, and
  \eqn{\gamma}-diversity.  \emph{Am. Nat.}, \bold{162}, 734--743.

  Lande, R.  (1996). Statistics and partitioning of species diversity,
  and similarity among multiple communities.  \emph{Oikos}, \bold{76},
  5--13.

}

\author{

  \enc{Péter Sólymos}{Peter Solymos}, \email{solymos@ualberta.ca}}
  \seealso{See \code{\link{oecosimu}} for permutation settings and
  calculating \eqn{p}-values.  }

\examples{
## NOTE: 'nsimul' argument usually needs to be >= 99
## here much lower value is used for demonstration

data(mite)
data(mite.xy)
data(mite.env)
## Function to get equal area partitions of the mite data
cutter <- function (x, cut = seq(0, 10, by = 2.5)) {
    out <- rep(1, length(x))
    for (i in 2:(length(cut) - 1))
        out[which(x > cut[i] & x <= cut[(i + 1)])] <- i
    return(out)}
## The hierarchy of sample aggregation
levsm <- data.frame(
    l1=1:nrow(mite),
    l2=cutter(mite.xy$y, cut = seq(0, 10, by = 2.5)),
    l3=cutter(mite.xy$y, cut = seq(0, 10, by = 5)),
    l4=cutter(mite.xy$y, cut = seq(0, 10, by = 10)))
## Let's see in a map
par(mfrow=c(1,3))
plot(mite.xy, main="l1", col=as.numeric(levsm$l1)+1)
plot(mite.xy, main="l2", col=as.numeric(levsm$l2)+1)
plot(mite.xy, main="l3", col=as.numeric(levsm$l3)+1)
par(mfrow=c(1,1))
## Additive diversity partitioning
adipart(mite, index="richness", nsimul=19)
adipart(mite ~ ., levsm, index="richness", nsimul=19)
## Hierarchical null model testing
## diversity analysis (similar to adipart)
hiersimu(mite, FUN=diversity, relative=TRUE, nsimul=19)
hiersimu(mite ~., levsm, FUN=diversity, relative=TRUE, nsimul=19)
## Hierarchical testing with the Morisita index
morfun <- function(x) dispindmorisita(x)$imst
hiersimu(mite ~., levsm, morfun, drop.highest=TRUE, nsimul=19)
}
\keyword{multivariate}

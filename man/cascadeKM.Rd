\name{cascadeKM}
\alias{cascadeKM}
\alias{cIndexKM}
\alias{plot.cascadeKM}
\alias{orderingKM}
\alias{pregraphKM}

\title{K-means partitioning using a range of values of K}
\description{
 This function is a wrapper for the \code{kmeans} function. It creates
 several partitions forming a cascade from a small to a large number of
 groups. 
}

\usage{
cascadeKM(data, inf.gr, sup.gr, iter = 100, criterion = "calinski")

cIndexKM(y, x, index = "all")

\method{plot}{cascadeKM}(x, min.g, max.g, grpmts.plot = TRUE, 
     sortg = FALSE, gridcol = NA, ...) 
}

\arguments{
  \item{ data }{ The data matrix. The objects (samples) are the rows.}
  \item{ inf.gr }{ The number of groups for the partition with the 
    smallest number of groups of the cascade (min).}
  \item{ sup.gr }{ The number of groups for the partition with the largest 	
    number of groups of the cascade (max).}
  \item{ iter }{ The number of random starting configurations for each value
    of \eqn{K}.}
  \item{ criterion }{ The criterion that will be used to select the best
    partition. The default value is \code{"calinski"}, which refers to
    the Calinski-Harabasz (1974) criterion. The simple structure index
    (\code{"ssi"}) is also available. Other indices are available in
    function \code{\link[cclust]{clustIndex}} (package \pkg{cclust}). In
    our experience, the two indices that work best and are most likely
    to return their maximum value at or near the optimal number of
    clusters are \code{"calinski"} and \code{"ssi"}. }
  \item{y}{Object of class \code{"kmeans"} returned by a clustering algorithm
    such as \code{\link{kmeans}}}
  \item{x}{Data matrix where columns correspond to variables and rows to
    observations, or the plotting object in \code{plot}}
  \item{index}{The available indices are: \code{"calinski"} and \code{"ssi"}. 
    Type \code{"all"} to obtain both indices. 
    Abbreviations of these names are also accepted.}
  \item{min.g, max.g}{The minimum and maximum numbers of groups to be
    displayed.}
  \item{grpmts.plot}{Show the plot (\code{TRUE} or \code{FALSE}).}
  \item{sortg}{Sort the objects as a function of their group membership
    to produce a more easily interpretable graph. See Details. The
    original object names are kept; they are used as labels in the
    output table \code{x}, although not in the graph.  If there were no
    row names, sequential row numbers are used to keep track of the
    original order of the objects.}
  \item{gridcol}{The colour of the grid lines in the plots. \code{NA},
    which is the default value, removes the grid lines.}
  \item{\dots}{Other parameters to the functions (ignored).}
    
}
\details{
  The function creates several partitions forming a cascade from a small
  to a large number of groups formed by \code{\link{kmeans}}.  Most
  of the work is performed by function \code{cIndex} which is based on the
  \code{\link[cclust]{clustIndex}} function (package \pkg{cclust}). 
  Some of the criteria were removed from this version because computation 
  errors were generated when only one object was found in a group.
  
  The default value is \code{"calinski"}, which refers to the well-known
  Calinski-Harabasz (1974) criterion. The other available index is the
  simple structure index \code{"ssi"} (Dolnicar et al. 1999).
  In the case of groups of equal
  sizes, \code{"calinski"} is generally a good criterion to indicate the
  correct number of groups. Users should not take its indications
  literally when the groups are not equal in size. Type \code{"all"} to
  obtain  both indices. The indices are defined as: 
  \describe{
    \item{calinski:}{
    \eqn{(SSB/(K-1))/(SSW/(n-K))}, where \eqn{n} is the
    number of data points and \eqn{K} is the number of clusters.
    \eqn{SSW} is the sum of squares within the clusters while
    \eqn{SSB} is the sum of squares among the clusters. This index
    is simply an \eqn{F} (ANOVA) statistic.}
 
    \item{ssi:}{
    the \dQuote{Simple Structure Index} multiplicatively combines
    several elements which influence the interpretability of a
    partitioning solution. The best partition is indicated by the
    highest SSI value.}
  }

  In a simulation study, Milligan and Cooper (1985) found
  that the Calinski-Harabasz criterion recovered the correct number of
  groups the most often. We recommend this criterion because, if the
  groups are of equal sizes, the maximum value of \code{"calinski"}
  usually indicates the correct number of groups. Another available
  index is the simple structure index \code{"ssi"}. Users should not
  take the indications of these indices literally when the groups are
  not equal in size and explore the groups corresponding to other values
  of \eqn{K}.
  
  Function \code{cascadeKM} has a \code{plot} method.  Two plots are
  produced. The graph on the left has the objects in 
  abscissa and the number of groups in ordinate. The groups are
  represented by colours. The graph on the right shows the values of the
  criterion (\code{"calinski"} or \code{"ssi"}) for determining the best
  partition. The highest value of the criterion is marked in red. Points
  marked in orange, if any, indicate partitions producing an increase in
  the criterion value as the number of groups increases; they may
  represent other interesting partitions.
  
  If \code{sortg=TRUE}, the objects are reordered by the following
  procedure: (1) a simple matching distance matrix is computed among the
  objects, based on the table of K-means assignments to groups, from
  \eqn{K} = \code{min.g} to \eqn{K} = \code{max.g}. (2) A principal
  coordinate analysis (PCoA, Gower 1966) is computed on the centred
  distance matrix. (3) The first principal coordinate is used as the new
  order of the objects in the graph. A simplified algorithm is used to
  compute the first principal coordinate only, using the iterative
  algorithm described in Legendre & Legendre (2012). The
  full distance matrix among objects is never computed; this avoids
  the problem of storing it when the number of objects is
  large. Distance values are computed as they are needed by the
  algorithm.
}

\value{ Function \code{cascadeKM} returns an object of class
  \code{cascadeKM} with items:
  \item{ partition }{ Table with the partitions found for different numbers 
    of groups \eqn{K}, from \eqn{K} = \code{inf.gr} to \eqn{K} =
    \code{sup.gr}. } 
  \item{ results }{ Values of the criterion to select the best
    partition. } 
  \item{ criterion }{ The name of the criterion used. }
  \item{ size }{ The number of objects found in each group, for all 
    partitions (columns). }

  Function \code{cIndex} returns a vector with the index values. The
  maximum value of these indices is supposed to indicate the best
  partition. These indices work best with groups of equal sizes. When
  the groups are not of equal sizes, one should not put too much faith
  in the maximum of these indices, and also explore the groups
  corresponding to other values of \eqn{K}.
}
\references{

  Calinski, T. and J. Harabasz. 1974. A dendrite method for cluster
  analysis. \emph{Commun. Stat.} \strong{3}: 1--27.

  Dolnicar, S., K. Grabler and J. A. Mazanec. 1999.  A tale of three
  cities: perceptual charting for analyzing destination images. Pp.
  39-62 in: Woodside, A. et al. [eds.] \emph{Consumer psychology of
  tourism, hospitality and leisure}. CAB International, New York.

  
  Gower, J. C. 1966. Some distance properties of latent root and vector
  methods used in multivariate analysis. \emph{Biometrika} \strong{53}:
  325--338.
  
  Legendre, P. & L. Legendre. 2012. \emph{Numerical ecology}, 3rd
  English edition. Elsevier Science BV, Amsterdam.
  
  Milligan, G. W. & M. C. Cooper. 1985. An examination of procedures for
  determining the number of clusters in a data set. \emph{Psychometrika}
  \strong{50}: 159--179.

  Weingessel, A., Dimitriadou, A. and Dolnicar, S. 2002. An examination
  of indexes for determining the number of clusters in binary data
  sets. \emph{Psychometrika} \strong{67}: 137--160.
}

\author{ Marie-Helene Ouellette
  \email{Marie-Helene.Ouellette@UMontreal.ca}, Sebastien Durand
  \email{Sebastien.Durand@UMontreal.ca} and Pierre Legendre
  \email{Pierre.Legendre@UMontreal.ca}. Edited for \pkg{vegan} by Jari
  Oksanen.
}

\seealso{\code{\link{kmeans}}, \code{\link[cclust]{clustIndex}}.}

\examples{
 # Partitioning a (10 x 10) data matrix of random numbers
 mat <- matrix(runif(100),10,10)
 res <- cascadeKM(mat, 2, 5, iter = 25, criterion = 'calinski') 
 toto <- plot(res)
 
 # Partitioning an autocorrelated time series
 vec <- sort(matrix(runif(30),30,1))
 res <- cascadeKM(vec, 2, 5, iter = 25, criterion = 'calinski')
 toto <- plot(res)
 
 # Partitioning a large autocorrelated time series
 # Note that we remove the grid lines
 vec <- sort(matrix(runif(1000),1000,1))
 res <- cascadeKM(vec, 2, 7, iter = 10, criterion = 'calinski')
 toto <- plot(res, gridcol=NA)
 
}
\keyword{ cluster }

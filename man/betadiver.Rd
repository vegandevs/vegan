\encoding{UTF-8}
\name{betadiver}
\alias{betadiver}
\alias{scores.betadiver}
\alias{plot.betadiver}

\title{ Indices of beta Diversity }
\description{
  The function estimates any of the 24 indices of beta diversity
  reviewed by Koleff et al. (2003). Alternatively, it finds the
  co-occurrence frequencies for triangular plots (Koleff et
  al. 2003). 
}

\usage{
betadiver(x, method = NA, order = FALSE, help = FALSE, ...)
\method{plot}{betadiver}(x, ...)
\method{scores}{betadiver}(x, triangular = TRUE, ...)
}

\arguments{
  \item{x}{Community data matrix, or the \code{betadiver} result for
  \code{plot} and \code{scores} functions. }
  \item{method}{The index of beta diversity as defined in Koleff et al.
  (2003), Table 1. You can use either the subscript of \eqn{\beta} or
  the number of the index. See argument \code{help} below. }
  \item{order}{Order sites by increasing number of species. This will
  influence the configuration in the triangular plot and non-symmetric
  indices. }
  \item{help}{Show the numbers, subscript names and the defining
    equations of the indices and exit.}
  \item{triangular}{Return scores suitable for triangular plotting of
  proportions. If \code{FALSE}, returns a 3-column matrix of raw counts.}
  \item{\dots}{ Other arguments to functions. }
}

\details{ The most commonly used index of beta diversity is
  \eqn{\beta_w = S/\alpha - 1}, where \eqn{S} is the total number of
  species, and \eqn{\alpha} is the average number of species per site
  (Whittaker 1960). A drawback of this model is that \eqn{S} increases
  with sample size, but the expectation of \eqn{\alpha} remains
  constant, and so the beta diversity increases with sample size. A
  solution to this problem is to study the beta diversity of pairs of
  sites (Marion et al. 2017). If we denote the number of species
  shared between two sites as \eqn{a} and the numbers of unique
  species (not shared) as \eqn{b} and \eqn{c}, then \eqn{S = a + b +
  c} and \eqn{\alpha = (2 a + b + c)/2} so that \eqn{\beta_w =
  (b+c)/(2 a + b + c)}. This is the \enc{Sørensen}{Sorensen}
  dissimilarity as defined in \pkg{vegan} function
  \code{\link{vegdist}} with argument \code{binary = TRUE}. Many other
  indices are dissimilarity indices as well.

  Function \code{betadiver} finds all indices reviewed by Koleff et
  al. (2003). All these indices could be found with function
  \code{\link{designdist}}, but the current function provides a
  conventional shortcut. The function only finds the indices. The proper
  analysis must be done with functions such as \code{\link{betadisper}},
  \code{\link{adonis2}} or \code{\link{mantel}}.

  The indices are directly taken from Table 1 of Koleff et al. (2003),
  and they can be selected either by the index number or the subscript
  name used by Koleff et al. The numbers, names and defining equations
  can be seen using \code{betadiver(help = TRUE)}. In all cases where
  there are two alternative forms, the one with the term \eqn{-1} is
  used. There are several duplicate indices, and the number of distinct
  alternatives is much lower than 24 formally provided. The formulations
  used in functions differ occasionally from those in Koleff et
  al. (2003), but they are still mathematically equivalent. With
  \code{method = NA}, no index is calculated, but instead an object of
  class \code{betadiver} is returned. This is a list of elements
  \code{a}, \code{b} and \code{c}. Function \code{plot} can be used to
  display the proportions of these elements in triangular plot as
  suggested by Koleff et al. (2003), and \code{scores} extracts the
  triangular coordinates or the raw scores. Function \code{plot} returns
  invisibly the triangular coordinates as an \code{"\link{ordiplot}"}
  object. 
}

\value{ With \code{method = NA}, the function returns an object of class
  \code{"betadisper"} with elements \code{a}, \code{b}, and \code{c}. If
  \code{method} is specified, the function returns a \code{"dist"}
  object which can be used in any function analysing
  dissimilarities. For beta diversity, particularly useful functions are
  \code{\link{betadisper}} to study the betadiversity in groups,
  \code{\link{adonis2}} for any model, and \code{\link{mantel}} to
  compare beta diversities to other dissimilarities or distances
  (including geographical distances). Although \code{betadiver} returns
  a \code{"dist"} object, some indices are similarities and cannot be
  used as such in place of dissimilarities, but that is a user
  error. Functions 10 (\code{"j"}), 11 (\code{"sor"}) and 21
  (\code{"rlb"}) are similarity indices. Function sets argument
  \code{"maxdist"} similarly as \code{\link{vegdist}}, using \code{NA}
  when there is no fixed upper limit, and 0 for similarities.
}

\references{

  Baselga, A. (2010) Partitioning the turnover and nestedness
  components of beta diversity. \emph{Global Ecology and Biogeography}
  19, 134--143.

  Koleff, P., Gaston, K.J. and Lennon, J.J. (2003) Measuring beta
  diversity for presence-absence data. \emph{Journal of Animal
  Ecology} 72, 367--382.

  Marion, Z.H., Fordyce, J.A. and Fitzpatrick, B.M. (2017) Pairwise
  beta diversity resolves an underappreciated source of confusion in
  calculating species turnover. \emph{Ecology} 98, 933--939.

  Whittaker, R.H. (1960) Vegetation of Siskiyou mountains, Oregon and
  California. \emph{Ecological Monographs} 30, 279--338.

 }

\author{Jari Oksanen }
\section{Warning }{Some indices return similarities instead of dissimilarities.} 

\seealso{ \code{\link{designdist}} can be used to implement all these
  functions, and also allows using notation with \code{alpha} and
  \code{gamma} diversities.  \code{\link{vegdist}} has some canned
  alternatives.  Functions \code{\link{betadisper}},
  \code{\link{adonis2}} and \code{\link{mantel}} can be used for
  analysing beta diversity objects. The returned dissimilarities can
  be used in any distance-based methods, such as
  \code{\link{metaMDS}}, \code{\link{capscale}} and
  \code{\link{dbrda}}. Functions \code{\link{nestedbetasor}} and
  \code{\link{nestedbetajac}} implement decomposition beta diversity
  measures (\enc{Sørensen}{Sorensen} and Jaccard) into turnover and
  nestedness components following Baselga (2010).  }


\examples{
## Raw data and plotting
data(sipoo)
m <- betadiver(sipoo)
plot(m)
## The indices
betadiver(help=TRUE)
## The basic Whittaker index
d <- betadiver(sipoo, "w")
## This should be equal to Sorensen index (binary Bray-Curtis in
## vegan)
range(d - vegdist(sipoo, binary=TRUE))
}

\keyword{ multivariate }


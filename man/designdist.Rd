\encoding{UTF-8}
\name{designdist}
\alias{designdist}

\title{Design your own Dissimilarities }
\description{
  You can define your own dissimilarities using terms for shared and
  total quantities, number of rows and number of columns. The shared and
  total quantities can be binary, quadratic or minimum terms. In
  binary terms, the shared component is number of shared species, and
  totals are numbers of species on sites. The quadratic terms
  are cross-products and sums of squares, and minimum terms 
  are sums of parallel minima and row totals.
}
\usage{
designdist(x, method = "(A+B-2*J)/(A+B)",
           terms = c("binary", "quadratic", "minimum"), 
           abcd = FALSE, name)
}

\arguments{
  \item{x}{Input data. }
  \item{method}{Equation for your dissimilarities. This can use terms
    \code{J} for shared quantity, \code{A} and \code{B} for totals,
    \code{N} for the number of rows (sites) and \code{P} for the number of
    columns (species). The equation can also contain any \R functions
    that accepts vector arguments and returns vectors of the same
    length. }
  \item{terms}{How shared and total components are found. For vectors
    \code{x} and \code{y} the  \code{"quadratic"} terms are \code{J = sum(x*y)},
    \code{A = sum(x^2)}, \code{B = sum(y^2)}, and \code{"minimum"} terms
    are \code{J = sum(pmin(x,y))}, \code{A = sum(x)} and \code{B =
      sum(y)}, and \code{"binary"} terms are either of these after transforming
    data into binary form (shared number of species, and number of
    species for each row). }
  \item{abcd}{Use 2x2 contingency table notation for binary data:
  \eqn{a} is the number of shared species, \eqn{b} and \eqn{c} are the
  numbers of species occurring only one of the sites but not in both,
  and \eqn{d} is the number of species that occur on neither of the sites.}
  \item{name}{The name you want to use for your index. The default is to
    combine the \code{method} equation and \code{terms} argument.}
}
\details{
  Most popular dissimilarity measures in ecology can be expressed with
  the help of terms \code{J}, \code{A} and \code{B}, and some also involve
  matrix dimensions \code{N} and \code{P}. Some examples you can define in
  \code{designdist} are:
  \tabular{lll}{
    \code{A+B-2*J} \tab \code{"quadratic"} \tab squared Euclidean \cr
    \code{A+B-2*J} \tab \code{"minimum"} \tab Manhattan \cr
    \code{(A+B-2*J)/(A+B)} \tab \code{"minimum"} \tab Bray-Curtis \cr
    \code{(A+B-2*J)/(A+B)} \tab \code{"binary"} \tab
    \enc{Sørensen}{Sorensen} \cr 
    \code{(A+B-2*J)/(A+B-J)} \tab \code{"binary"} \tab Jaccard \cr
    \code{(A+B-2*J)/(A+B-J)} \tab \code{"minimum"} \tab
    \enc{Ružička}{Ruzicka} \cr
    \code{(A+B-2*J)/(A+B-J)} \tab \code{"quadratic"} \tab
    (dis)similarity ratio \cr
    \code{1-J/sqrt(A*B)} \tab \code{"binary"} \tab Ochiai \cr
    \code{1-J/sqrt(A*B)} \tab \code{"quadratic"} \tab cosine
    complement \cr
    \code{1-phyper(J-1, A, P-A, B)} \tab \code{"binary"} \tab Raup-Crick 
  }

  The function \code{designdist} can implement most dissimilarity
  indices in \code{\link{vegdist}} or elsewhere, and it can also be
  used to implement many other indices, amongst them, most of those
  described in Legendre & Legendre (1998). It can also be used to
  implement all indices of beta diversity described in Koleff et
  al. (2003), but there also is a specific function
  \code{\link{betadiver}} for the purpose.

  If you want to implement binary dissimilarities based on the 2x2
  contingency table notation, you can set \code{abcd = TRUE}. In this
  notation \code{a = J}, \code{b = A-J}, \code{c = B-J}, \code{d =
  P-A-B+J}. This notation is often used instead fo the more more
  tangible default notation for reasons that are opaque to me. 
}

\value{
  \code{designdist} returns an object of class \code{\link{dist}}.
}
\references{
  Koleff, P., Gaston, K.J. and Lennon, J.J. (2003) Measuring beta
  diversity for presence--absence data. \emph{J. Animal Ecol.}
  \strong{72}, 367--382. 
  
  Legendre, P. and Legendre, L. (1998) \emph{Numerical Ecology}. 2nd
  English ed. Elsevier
  }
\author{ Jari Oksanen }
\note{  
  \code{designdist} does not use compiled code, and may be slow or use
  plenty of memory in large data sets. It is very easy to make errors
  when  defining a function by hand. If an index is available in a
  function  using compiled code, it is better to use the canned
  alternative.
}
\seealso{ \code{\link{vegdist}}, \code{\link{betadiver}}, \code{\link{dist}}. }
\examples{
## Arrhenius dissimilarity: the value of z in the species-area model
## S = c*A^z when combining two sites of equal areas, where S is the
## number of species, A is the area, and c and z are model parameters.
## The A below is not the area (which cancels out), but number of
## species in one of the sites, as defined in designdist().
data(BCI)
dis <- designdist(BCI, "(log(A+B-J)-log(A+B)+log(2))/log(2)")
## This can be used in clustering or ordination...
ordiplot(cmdscale(dis))
## ... or in analysing beta diversity (without gradients)
summary(dis)
  }

\keyword{multivariate }


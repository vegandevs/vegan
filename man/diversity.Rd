\encoding{UTF-8}
\name{diversity}
\alias{diversity}
\alias{fisher.alpha}
\alias{specnumber}

\title{Ecological Diversity Indices}

\description{ Shannon, Simpson, and Fisher diversity indices and species
  richness.  }

\usage{
diversity(x, index = "shannon", MARGIN = 1, base = exp(1))
fisher.alpha(x, MARGIN = 1, ...)
specnumber(x, groups, MARGIN = 1)
}

\arguments{
  \item{x}{Community data, a matrix-like object or a vector.}
  \item{index}{Diversity index, one of \code{"shannon"},
    \code{"simpson"} or \code{"invsimpson"}.}
  \item{MARGIN}{Margin for which the index is computed. }
  \item{base}{ The logarithm \code{base} used in \code{shannon}.}
  \item{groups}{A grouping factor: if given, finds the total number of
    species in each group.}
  \item{...}{Parameters passed to the function.}
}
\details{
  Shannon or Shannon--Weaver (or Shannon--Wiener) index is defined as
  \eqn{H' = -\sum_i p_i \log_{b} p_i}{H = -sum p_i log(b) p_i}, where
  \eqn{p_i} is the proportional abundance of species \eqn{i} and \eqn{b}
  is the base of the logarithm.  It is most popular to use natural
  logarithms, but some argue for base \eqn{b = 2} (which makes sense,
  but no real difference).

  Both variants of Simpson's index are based on \eqn{D = \sum p_i^2}{D =
    sum p_i^2}. Choice \code{simpson} returns \eqn{1-D} and
  \code{invsimpson} returns \eqn{1/D}.
 
  \code{fisher.alpha} estimates the \eqn{\alpha} parameter of
  Fisher's logarithmic series (see \code{\link{fisherfit}}). 
  The estimation is possible only for genuine
  counts of individuals. 
  
  Function \code{specnumber} finds the number of species. With
  \code{MARGIN = 2}, it finds frequencies of species. If \code{groups}
  is given, finds the total number of species in each group (see 
  example on finding one kind of beta diversity with this option).
  
  Better stories can be told about Simpson's index than about
  Shannon's index, and still grander narratives about
  rarefaction (Hurlbert 1971).  However, these indices are all very
  closely related (Hill 1973), and there is no reason to despise one
  more than others (but if you are a graduate student, don't drag me in,
  but obey your Professor's orders). In particular, the exponent of the
  Shannon index is linearly related to inverse Simpson (Hill 1973)
  although the former may be more sensitive to rare species. Moreover,
  inverse Simpson is asymptotically equal to rarefied species richness
  in sample of two individuals, and Fisher's \eqn{\alpha} is very
  similar to inverse Simpson.

}

\value{
  A vector of diversity indices or numbers of species. 
}

\references{
  Fisher, R.A., Corbet, A.S. & Williams, C.B. (1943). The relation
  between the number of species and the number of individuals in a
  random sample of animal population. \emph{Journal of Animal Ecology}
  \strong{12}, 42--58.
  
  Hurlbert, S.H. (1971). The nonconcept of species diversity: a critique
  and alternative parameters. \emph{Ecology} \strong{52}, 577--586.
 
}

\seealso{These functions calculate only some basic indices, but many
  others can be derived with them (see Examples). Facilities related to
  diversity are discussed in a \pkg{vegan} vignette that can be read
  with command \code{vegandocs("diversity")}.  Functions
  \code{\link{renyi}} and \code{\link{tsallis}} estimate a series of
  generalized diversity indices. Function \code{\link{rarefy}} finds
  estimated number of species for given sample size. Beta diversity can
  be estimated with \code{\link{betadiver}}. Diversities can be
  partitioned with \code{\link{adipart}} and \code{\link{multipart}}. }

\author{ Jari Oksanen  and  Bob O'Hara (\code{fisher.alpha}).}

\examples{
data(BCI)
H <- diversity(BCI)
simp <- diversity(BCI, "simpson")
invsimp <- diversity(BCI, "inv")
## Unbiased Simpson (Hurlbert 1971, eq. 5) with rarefy:
unbias.simp <- rarefy(BCI, 2) - 1
## Fisher alpha
alpha <- fisher.alpha(BCI)
## Plot all
pairs(cbind(H, simp, invsimp, unbias.simp, alpha), pch="+", col="blue")
## Species richness (S) and Pielou's evenness (J):
S <- specnumber(BCI) ## rowSums(BCI > 0) does the same...
J <- H/log(S)
## beta diversity defined as gamma/alpha - 1:
data(dune)
data(dune.env)
alpha <- with(dune.env, tapply(specnumber(dune), Management, mean))
gamma <- with(dune.env, specnumber(dune, Management))
gamma/alpha - 1
}
\keyword{ univar }



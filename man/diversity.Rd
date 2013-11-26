\encoding{UTF-8}
\name{diversity}
\alias{diversity}
\alias{rarefy}
\alias{rrarefy}
\alias{drarefy}
\alias{rarecurve}
\alias{fisher.alpha}
\alias{specnumber}

\title{ Ecological Diversity Indices and Rarefaction Species Richness }

\description{ Shannon, Simpson, and Fisher diversity indices
and rarefied species richness for community ecologists.  }

\usage{
diversity(x, index = "shannon", MARGIN = 1, base = exp(1))
rarefy(x, sample, se = FALSE, MARGIN = 1)
rrarefy(x, sample)
drarefy(x, sample)
rarecurve(x, step = 1, sample, xlab = "Sample Size", ylab = "Species",
   label = TRUE, ...)
fisher.alpha(x, MARGIN = 1, ...)
specnumber(x, groups, MARGIN = 1)
}

\arguments{
  \item{x}{Community data, a matrix-like object or a vector.}
  \item{index}{Diversity index, one of \code{"shannon"},
    \code{"simpson"} or \code{"invsimpson"}.}
  \item{MARGIN}{Margin for which the index is computed. }
  \item{base}{ The logarithm \code{base} used in \code{shannon}.}
  \item{sample}{Subsample size for rarefying community, either a single
    value or a vector.}
  \item{se}{Estimate standard errors.}
  \item{step}{Step size for sample sizes in rarefaction curves.}
  \item{xlab, ylab}{Axis labels in plots of rarefaction curves.}
  \item{label}{Label rarefaction curves by rownames of \code{x} (logical).}
  \item{groups}{A grouping factor: if given, finds the total number of
    species in each group.}
  \item{...}{Parameters passed to \code{\link{nlm}}, or to \code{\link{plot}}, 
    \code{\link{lines}} and \code{\link{ordilabel}} in \code{rarecurve}.}
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
 
  Function \code{rarefy} gives the expected species richness in random
  subsamples of size \code{sample} from the community. The size of
  \code{sample} should be smaller than total community size, but the
  function will silently work for larger \code{sample} as well and
  return non-rarefied species richness (and standard error = 0). If
  \code{sample} is a vector, rarefaction of all observations is
  performed for each sample size separately.  Rarefaction can be
  performed only with genuine counts of individuals.  The function
  \code{rarefy} is based on Hurlbert's (1971) formulation, and the
  standard errors on Heck et al. (1975).

  Function \code{rrarefy} generates one randomly rarefied community data
  frame or vector of given \code{sample} size. The \code{sample} can be
  a vector giving the sample sizes for each row, and its values must be
  less or equal to observed number of individuals. The random
  rarefaction is made without replacement so that the variance of
  rarefied communities is rather related to rarefaction proportion than
  to to the size of the \code{sample}.

  Function \code{drarefy} returns probabilities that species occur in a
  rarefied community of size \code{sample}. The \code{sample} can be a
  vector giving the sample sizes for each row.

  Function \code{rarecurve} draws a rarefaction curve for each row of
  the input data. The rarefaction curves are evaluated using the
  interval of \code{step} sample sizes, always including 1 and total
  sample size.  If \code{sample} is specified, a vertical line is
  drawn at \code{sample} with horizontal lines for the rarefied
  species richnesses.

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
  A vector of diversity indices or rarefied species richness values. With
  a single \code{sample} and \code{se = TRUE}, function \code{rarefy}
  returns a 2-row matrix 
  with rarefied richness (\code{S}) and its standard error
  (\code{se}). If \code{sample} is a vector in \code{rarefy}, the
  function returns a matrix with a column for each \code{sample} size,
  and if \code{se = TRUE}, rarefied richness and its standard error are
  on consecutive lines.

  Function \code{rarecurve} returns \code{\link{invisible}} list of
  \code{rarefy} results corresponding each drawn curve.

  With option \code{se = TRUE}, function \code{fisher.alpha} returns a
  data frame with items for \eqn{\alpha} (\code{alpha}), its approximate
  standard errors (\code{se}), residual degrees of freedom
  (\code{df.residual}), and the \code{code} returned by
  \code{\link{nlm}} on the success of estimation. 
}

\references{
  Fisher, R.A., Corbet, A.S. & Williams, C.B. (1943). The relation
  between the number of species and the number of individuals in a
  random sample of animal population. \emph{Journal of Animal Ecology}
  \strong{12}, 42--58.

  Heck, K.L., van Belle, G. & Simberloff, D. (1975). Explicit
  calculation of the rarefaction diversity measurement and the
  determination of sufficient sample size. \emph{Ecology} \strong{56},
  1459--1461.  
  
  Hurlbert, S.H. (1971). The nonconcept of species diversity: a critique
  and alternative parameters. \emph{Ecology} \strong{52}, 577--586.
 
}

\seealso{Function \code{renyi} for generalized \enc{Rényi}{Renyi}
  diversity and Hill numbers.}

\author{ Jari Oksanen  and  Bob O'Hara (\code{fisher.alpha}).}

\examples{
data(BCI)
H <- diversity(BCI)
simp <- diversity(BCI, "simpson")
invsimp <- diversity(BCI, "inv")
## Unbiased Simpson of Hurlbert 1971 (eq. 5):
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
## Rarefaction
(raremax <- min(rowSums(BCI)))
Srare <- rarefy(BCI, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(BCI, step = 20, sample = raremax, col = "blue", cex = 0.6)
}
\keyword{ univar }



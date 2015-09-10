\name{rarefy}
\alias{rarefy}
\alias{rrarefy}
\alias{drarefy}
\alias{rarecurve}
\alias{rareslope}

\title{Rarefaction Species Richness}

\description{ Rarefied species richness for community ecologists.  }

\usage{
rarefy(x, sample, se = FALSE, MARGIN = 1)
rrarefy(x, sample)
drarefy(x, sample)
rarecurve(x, step = 1, sample, xlab = "Sample Size", ylab = "Species",
          label = TRUE, col, lty, ...)
rareslope(x, sample)
}

\arguments{
  \item{x}{Community data, a matrix-like object or a vector.}
  \item{MARGIN}{Margin for which the index is computed. }
  \item{sample}{Subsample size for rarefying community, either a single
    value or a vector.}
  \item{se}{Estimate standard errors.}
  \item{step}{Step size for sample sizes in rarefaction curves.}
  \item{xlab, ylab}{Axis labels in plots of rarefaction curves.}
  \item{label}{Label rarefaction curves by rownames of \code{x}
    (logical).}
  \item{col, lty}{plotting colour and line type, see
    \code{\link{par}}. Can be a vector of length \code{nrow(x)}, one per
    sample, and will be extended to such a length internally.}
  \item{...}{Parameters passed to \code{\link{nlm}}, or to \code{\link{plot}}, 
    \code{\link{lines}} and \code{\link{ordilabel}} in \code{rarecurve}.}
}
\details{
 
  Function \code{rarefy} gives the expected species richness in random
  subsamples of size \code{sample} from the community. The size of
  \code{sample} should be smaller than total community size, but the
  function will work for larger \code{sample} as well (with a warning)
  and return non-rarefied species richness (and standard error =
  0). If \code{sample} is a vector, rarefaction of all observations is
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

  Function \code{rareslope} calculates the slope of \code{rarecurve}
  (derivative of \code{rarefy}) at given \code{sample} size; the
  \code{sample} need not be an integer.

}

\value{
  A vector of rarefied species richness values. With a single
  \code{sample} and \code{se = TRUE}, function \code{rarefy} returns a
  2-row matrix with rarefied richness (\code{S}) and its standard error
  (\code{se}). If \code{sample} is a vector in \code{rarefy}, the
  function returns a matrix with a column for each \code{sample} size,
  and if \code{se = TRUE}, rarefied richness and its standard error are
  on consecutive lines.

  Function \code{rarecurve} returns \code{\link{invisible}} list of
  \code{rarefy} results corresponding each drawn curve.
}

\references{
  Heck, K.L., van Belle, G. & Simberloff, D. (1975). Explicit
  calculation of the rarefaction diversity measurement and the
  determination of sufficient sample size. \emph{Ecology} \strong{56},
  1459--1461.  
  
  Hurlbert, S.H. (1971). The nonconcept of species diversity: a critique
  and alternative parameters. \emph{Ecology} \strong{52}, 577--586.
 
}

\seealso{Use \code{\link{specaccum}} for species accumulation curves
  where sites are sampled instead of individuals. \code{\link{specpool}}
  extrapolates richness to an unknown sample size.}

\author{Jari Oksanen}

\examples{
data(BCI)
S <- specnumber(BCI) # observed number of species
(raremax <- min(rowSums(BCI)))
Srare <- rarefy(BCI, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(BCI, step = 20, sample = raremax, col = "blue", cex = 0.6)
}
\keyword{ univar }



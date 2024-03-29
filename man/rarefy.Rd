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
          label = TRUE, col, lty, tidy = FALSE, ...)
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
  \item{tidy}{Instead of drawing a \code{plot}, return a \dQuote{tidy}
    data frame than can be used in \CRANpkg{ggplot2} graphics. The data
    frame has variables \code{Site} (factor), \code{Sample} and
    \code{Species}.}
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

  Function \code{rrarefy} generates one randomly rarefied community
  data frame or vector of given \code{sample} size. The \code{sample}
  can be a vector giving the sample sizes for each row.  If the
  \code{sample} size is equal to or larger than the observed number
  of individuals, the non-rarefied community will be returned.  The
  random rarefaction is made without replacement so that the variance
  of rarefied communities is rather related to rarefaction proportion
  than to the size of the \code{sample}. Random rarefaction is
  sometimes used to remove the effects of different sample
  sizes. This is usually a bad idea: random rarefaction discards valid
  data, introduces random error and reduces the quality of the data
  (McMurdie & Holmes 2014). It is better to use normalizing
  transformations (\code{\link{decostand}} in \pkg{vegan}) possible
  with variance stabilization (\code{\link{decostand}} and
  \code{\link{dispweight}} in \pkg{vegan}) and methods that are not
  sensitive to sample sizes.

  Function \code{drarefy} returns probabilities that species occur in
  a rarefied community of size \code{sample}. The \code{sample} can be
  a vector giving the sample sizes for each row. If the \code{sample}
  is equal to or larger than the observed number of individuals, all
  observed species will have sampling probability 1.

  Function \code{rarecurve} draws a rarefaction curve for each row of
  the input data. The rarefaction curves are evaluated using the
  interval of \code{step} sample sizes, always including 1 and total
  sample size.  If \code{sample} is specified, a vertical line is
  drawn at \code{sample} with horizontal lines for the rarefied
  species richnesses.

  Function \code{rareslope} calculates the slope of \code{rarecurve}
  (derivative of \code{rarefy}) at given \code{sample} size; the
  \code{sample} need not be an integer.

  Rarefaction functions should be used for observed counts. If you
  think it is necessary to use a multiplier to data, rarefy first and
  then multiply. Removing rare species before rarefaction can also
  give biased results. Observed count data normally include singletons
  (species with count 1), and if these are missing, functions issue
  warnings. These may be false positives, but it is recommended to
  check that the observed counts are not multiplied or rare taxa are
  not removed.

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
  \code{rarefy} results corresponding each drawn curve. Alternatively,
  with \code{tidy = TRUE} it returns a data frame that can be used in
  \CRANpkg{ggplot2} graphics.
}

\references{
  Heck, K.L., van Belle, G. & Simberloff, D. (1975). Explicit
  calculation of the rarefaction diversity measurement and the
  determination of sufficient sample size. \emph{Ecology} \strong{56},
  1459--1461.

  Hurlbert, S.H. (1971). The nonconcept of species diversity: a critique
  and alternative parameters. \emph{Ecology} \strong{52}, 577--586.

  McMurdie, P.J. & Holmes, S. (2014). Waste not, want not: Why
  rarefying microbiome data is inadmissible. \emph{PLoS Comput Biol}
  \strong{10(4):} e1003531. \doi{10.1371/journal.pcbi.1003531}

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



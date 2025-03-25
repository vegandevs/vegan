\encoding{UTF-8}
\name{dispindmorisita}
\alias{dispindmorisita}
\title{Morisita index of intraspecific aggregation}
\description{
Calculates the Morisita index of dispersion, standardized index values, and the so called clumpedness and uniform indices.
}
\usage{
dispindmorisita(x, unique.rm = FALSE, crit = 0.05, na.rm = FALSE)
}
\arguments{
  \item{x}{community data matrix, with sites (samples) as rows and
    species as columns.}
  \item{unique.rm}{logical, if \code{TRUE}, unique species (occurring
    in only one sample) are removed from the result.}
  \item{crit}{two-sided p-value used to calculate critical
    Chi-squared values.}
  \item{na.rm}{logical.
    Should missing values (including \code{NaN}) be omitted from the
    calculations?}
}
\details{
The Morisita index of dispersion is defined as (Morisita 1959, 1962):

\code{Imor = n * (sum(xi^2) - sum(xi)) / (sum(xi)^2 - sum(xi))}

where \eqn{xi} is the count of individuals in sample \eqn{i}, and
\eqn{n} is the number of samples (\eqn{i = 1, 2, \ldots, n}).
\eqn{Imor} has values from 0 to \eqn{n}. In uniform (hyperdispersed)
patterns its value falls between 0 and 1, in clumped patterns it falls
between 1 and \eqn{n}. For increasing sample sizes (i.e. joining
neighbouring quadrats), \eqn{Imor} goes to \eqn{n} as the
quadrat size approaches clump size. For random patterns,
\eqn{Imor = 1} and counts in the samples follow Poisson
frequency distribution.

The deviation from random expectation (null hypothesis)
can be tested using critical values of the Chi-squared
distribution with \eqn{n-1} degrees of freedom.
Confidence intervals around 1 can be calculated by the clumped
\eqn{Mclu} and uniform \eqn{Muni} indices (Hairston et al. 1971, Krebs
1999) (Chi2Lower and Chi2Upper refers to e.g. 0.025 and 0.975 quantile
values of the Chi-squared distribution with \eqn{n-1} degrees of
freedom, respectively, for \code{crit = 0.05}):

\code{Mclu = (Chi2Lower - n + sum(xi)) / (sum(xi) - 1)}

\code{Muni = (Chi2Upper - n + sum(xi)) / (sum(xi) - 1)}

Smith-Gill (1975) proposed scaling of Morisita index from [0, n]
interval into [-1, 1], and setting up -0.5 and 0.5 values as
confidence limits around random distribution with rescaled value 0. To
rescale the Morisita index, one of the following four equations apply
to calculate the standardized index \eqn{Imst}:

(a) \code{Imor >= Mclu > 1}: \code{Imst = 0.5 + 0.5 (Imor - Mclu) / (n - Mclu)},

(b) \code{Mclu > Imor >= 1}: \code{Imst = 0.5 (Imor - 1) / (Mclu - 1)},

(c) \code{1 > Imor > Muni}: \code{Imst = -0.5 (Imor - 1) / (Muni - 1)},

(d) \code{1 > Muni > Imor}: \code{Imst = -0.5 + 0.5 (Imor - Muni) / Muni}.
}

\value{ Returns a data frame with as many rows as the number of columns
in the input data, and with four columns. Columns are: \code{imor} the
unstandardized Morisita index, \code{mclu} the clumpedness index,
\code{muni} the uniform index, \code{imst} the standardized Morisita
index, \code{pchisq} the Chi-squared based probability for the null
hypothesis of random expectation.
}

\references{

Morisita, M. 1959. Measuring of the dispersion of individuals and
analysis of the distributional patterns.  \emph{Mem. Fac. Sci. Kyushu
Univ. Ser. E} 2, 215--235.

Morisita, M. 1962. Id-index, a measure of dispersion of individuals.
\emph{Res. Popul. Ecol.} 4, 1--7.

Smith-Gill, S. J. 1975. Cytophysiological basis of disruptive pigmentary
patterns in the leopard frog, \emph{Rana pipiens}. II.  Wild type and
mutant cell specific patterns. \emph{J. Morphol.} 146, 35--54.

Hairston, N. G., Hill, R. and Ritte, U. 1971. The interpretation of
aggregation patterns. In: Patil, G. P., Pileou, E. C. and Waters,
W. E. eds. \emph{Statistical Ecology 1: Spatial Patterns and Statistical
Distributions}. Penn. State Univ. Press, University Park.

Krebs, C. J. 1999. \emph{Ecological Methodology}. 2nd ed. Benjamin
Cummings Publishers.  }

\author{\enc{Péter Sólymos}{Peter Solymos}, \email{solymos@ualberta.ca}}

\note{ A common error found in several papers is that when
standardizing as in the case (b), the denominator is given as
\code{Muni - 1}. This results in a hiatus in the [0, 0.5] interval of
the standardized index. The root of this typo is the book of Krebs
(1999), see the Errata for the book (Page 217, currently
https://www.zoology.ubc.ca/~krebs/downloads/errors_2nd_printing.pdf).
}

\examples{
data(dune)
x <- dispindmorisita(dune)
x
y <- dispindmorisita(dune, unique.rm = TRUE)
y
dim(x) ## with unique species
dim(y) ## unique species removed
}
\keyword{ multivariate }
\keyword{ spatial }

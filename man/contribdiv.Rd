\encoding{UTF-8}
\name{contribdiv}
\alias{contribdiv}
\alias{plot.contribdiv}
\title{Contribution Diversity Approach}
\description{
The contribution diversity approach is based in the differentiation of within-unit and among-unit diversity by using additive diversity partitioning and unit distinctiveness.
}
\usage{
contribdiv(comm, index = c("richness", "simpson"),
     relative = FALSE, scaled = TRUE, drop.zero = FALSE)
\method{plot}{contribdiv}(x, sub, xlab, ylab, ylim, col, ...)
}
\arguments{
  \item{comm}{The community data matrix with samples as rows and species as column.}
  \item{index}{Character, the diversity index to be calculated.}
  \item{relative}{Logical, if \code{TRUE} then contribution diversity values are expressed as their signed deviation from their mean. See details.}
  \item{scaled}{Logical, if \code{TRUE} then relative contribution diversity values are scaled by the sum of gamma values (if \code{index = "richness"}) or by sum of gamma values times the number of rows in \code{comm} (if \code{index = "simpson"}). See details.}
  \item{drop.zero}{Logical, should empty rows dropped from the result? If empty rows are not dropped, their corresponding results will be \code{NA}s.}
  \item{x}{An object of class \code{"contribdiv"}.}
  \item{sub, xlab, ylab, ylim, col}{Graphical arguments passed to plot.}
  \item{\dots}{Other arguments passed to plot.}
}
\details{
This approach was proposed by Lu et al. (2007). Additive diversity partitioning (see \code{\link{adipart}} for more references) deals with the relation of mean alpha and the total (gamma) diversity. Although alpha diversity values often vary considerably. Thus, contributions of the sites to the total diversity are uneven. This site specific contribution is measured by contribution diversity components. A unit that has e.g. many unique species will contribute more to the higher level (gamma) diversity than another unit with the same number of species, but all of which common.

Distinctiveness of species \eqn{j} can be defined as the number of sites where it occurs (\eqn{n_j}), or the sum of its relative frequencies (\eqn{p_j}). Relative frequencies are computed sitewise and \eqn{sum_j{p_ij}}s at site \eqn{i} sum up to \eqn{1}.

The contribution of site \eqn{i} to the total diversity is given by \eqn{alpha_i = sum_j(1 / n_ij)} when dealing with richness and \eqn{alpha_i = sum(p_{ij} * (1 - p_{ij}))} for the Simpson index.

The unit distinctiveness of site \eqn{i} is the average of the species distinctiveness, averaging only those species which occur at site \eqn{i}. For species richness: \eqn{alpha_i = mean(n_i)} (in the paper, the second equation contains a typo, \eqn{n} is without index). For the Simpson index: \eqn{alpha_i = mean(n_i)}.

The Lu et al. (2007) gives an in-depth description of the different indices.
}
\value{
An object of class \code{"contribdiv"}.
}
\references{
Lu, H. P., Wagner, H. H. and Chen, X. Y. 2007. A contribution diversity approach to evaluate species diversity.
\emph{Basic and Applied Ecology}, 8, 1--12.
}
\author{\enc{P\'eter S\'olymos}{Peter Solymos}, \email{solymos@ualberta.ca}}
\seealso{\code{\link{adipart}}, \code{\link{diversity}}}
\examples{
## Artificial example given in
## Table 2 in Lu et al. 2007
x <- matrix(c(
1/3,1/3,1/3,0,0,0,
0,0,1/3,1/3,1/3,0,
0,0,0,1/3,1/3,1/3),
3, 6, byrow = TRUE,
dimnames = list(LETTERS[1:3],letters[1:6]))
x
## Compare results with Table 2
contribdiv(x, "richness")
contribdiv(x, "simpson")
contribdiv(x, "richness", relative = TRUE, scaled = FALSE)
contribdiv(x, "simpson", relative = TRUE, scaled = FALSE)
## BCI data set
data(BCI)
opar <- par(mfrow=c(2,2))
plot(contribdiv(BCI, "richness"), main = "Absolute")
plot(contribdiv(BCI, "richness", relative = TRUE), main = "Relative")
plot(contribdiv(BCI, "simpson"))
plot(contribdiv(BCI, "simpson", relative = TRUE))
par(opar)
}
\keyword{multivariate}

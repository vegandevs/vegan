\name{simulSADcomm}
\alias{simulSADcomm}
\title{
Simulate community matrix with constant SAD
}
\description{
This function simulates community matrices with the same species abundance distribution following patterns defined by a set of explanatory variables. This function was used to simulate community matrices in Blanchet et al. (2014)}
\usage{
simulSADcomm(sp.abund, expl.var, expl.rand.sel = TRUE, nexpl.comb = 2, 
             binary = FALSE, fix.expl = NULL, nsite = 50, weight = NULL, 
             range.weight = c(0, 2), sd.expl = FALSE, norm = c(0, 1))
}
\arguments{
  \item{sp.abund}{
A vector defining the number of species in a bin. See Details for more information.
}
  \item{expl.var}{
A matrix of explanatory variables to use to construct the species.
}
  \item{expl.rand.sel}{
Logical. Whether an explanatory variable should be randomly selected to construct a species (TRUE) or a fixed combination should be given (FALSE). (Default is TRUE)
}
  \item{nexpl.comb}{
Numeric. The number of explanatory variables that will be combined together to construct a species. Default is 2.
}
  \item{binary}{
Logical. Whether the site-by-species matrix presents abundances (FALSE) or presence/absence (TRUE). Default is FALSE.
}
  \item{fix.expl}{
A matrix that defines which combination of explanatory variables should be used to construct species. This argument is only active when expl.rand.sel=FALSE. See Details for more information.
}
  \item{nsite}{
Numeric. Number of sites (rows) in the resulting community matrix. See Details.
}
  \item{weight}{
A vector of regression coefficient used to give weight to each species. If NULL, weights are selected by randomly samping of a uniform distribution with a range defined by range.weight. Default is NULL.
}
  \item{range.weight}{
A vector of length 2 giving the minimum and maximum of a uniform distribution from which the \code{weight} will be sampled. Default is 0 and 2.
}
  \item{sd.expl}{
Logical. This argument is only active when \code{expl.rand.sel} is FALSE (That is when a fixed combination of explanatory variable is used to construct a community matrix). This could be either the standard deviation of the Normal error added when constructing a species is a multiplier of the standard deviation of the deterministic portion of the newly created explanatory variable (TRUE) or the pure standard deviation (FALSE). Default is FALSE.
}
  \item{norm}{
Vector of length 2 giving the mean and a multiplier of the standard deviation of the deterministic portion of a newly created explanatory variable. Default is mean = 0 and multiplier of the standard deviation of the new deterministic explanatory variable = 1.
}
}
\details{
The argument \code{sp.abund} defines the species-abundance distribution structure of the data following the binning proposed by Gray et al. (2006). For example, if the vector is (40,20,30), it means that there will be 40 species with 1 individual, 20 with 2 or 3 individuals, and 30 with 4 to 7 individuals.

The individuals are assigned to the sites according to the set of exlanatory variables given in \code{expl.var}. It is possible that a site occur with 0 individuals. They will be included in the community matrix and should be dealt with \emph{a posteriori}. 

If \code{expl.rand.sel} is TRUE, the explanatory variables are randomly sampled (without replacement) when combining (adding) explanatory variables together. The number of explanatory variables must be a multiple of nexpl.comb.

Error is included to a species by multiplying a weight to the explanatory variable used to construct the species and by adding a normally distributed error term to the same explanatory variable. An error term with a standard deviation equal to the standard deviation of the explanatory variable allows for the explanatory variable to explain roughly 50% of the constructed species.

\code{fix.expl} is a matrix that has as many rows as there are species and as many columns as \code{nexpl.comb} (number of explanatory variables to combine). The numbers in \code{fix.expl} are integers that refers to the columns of \code{expl.var}. When \code{fix.expl} is TRUE, \code{nexpl.comb} becomes meaningless.

If a presence-absence matrix is constructed (\code{binary}=TRUE), \code{sp.abund} should be constructed in such a way that no bin should include species with an abundance larger than the number of sites, otherwise an error message will be sent. Within, this constraint, if the maximum of the last bin (the one with the largest abundance) is larger than the number of site, it will be automatically changed to the number of sites-1.

This function was designed to do much more than the simulations generated in the work of Blanchet et al. (in press). It is meant to be used for future simulation studies.
}
\value{
\code{site.sp} : The site (rows) by species (column) community matrix generated.
\code{sel.expl} : A vector presenting the order explanatory variables used to construct each species. The order follows the order of the species.
}
\references{
Gray, J.S., A. Bjorgeaeter, and K.I. Ugland. 2006. On plotting species abundance distributions, \emph{Journal of Animal Ecology} \strong{75}:752--756.

Blanchet, F.G., P. Legendre, J.A.C. Bergeron, and F. He. 2014. Consensus RDA across dissimilarity coefficients for canonical ordination of community composition data. \emph{Ecological Monographs} \strong{84}: 491--511.
}
\author{
F. Guillaume Blanchet
}

\examples{
SAD<-c(1,2,4,6,4,2,1,0,0,0)
expl<-matrix(rnorm(400),ncol=8)
simulSADcomm(SAD,expl)
}
\keyword{ datagen }

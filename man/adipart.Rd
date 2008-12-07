\encoding{UTF-8}
\name{adipart}
\alias{adipart}
\alias{print.adipart}
\title{Additive Diversity Partitioning}
\description{
In additive diversity partitioning, mean values of alpha diversity at lower levels of a sampling hierarchy are compared to the total number of species in the entire data set (gamma diversity) in the form: gamma = mean(alpha) + beta. Thus beta = gamma - mean(alpha).
}
\usage{
adipart(matr, strata, index=c("richness", "shannon", "simpson"),
weights=c("unif", "prop"), relative = FALSE, nsimul=99, control, ...)
\method{print}{adipart}(x, ...)
}
\arguments{
  \item{matr}{A community data matrix with samples as rows and species as column.}
  \item{strata}{A matrix or data frame containing levels of sampling hierarchy as columns, columns from right to left will be treated as nested (first column is the lowest, last is the highest level).}
  \item{index}{Character, the diversity index to be calculated, one of  \code{c("richness", "shannon", "simpson")}.}
  \item{weights}{Character, \code{"unif"} for uniform weights, \code{"prop"} for weighting proportional to sample abundances to use in weighted averaging of individual alpha values within strata of a given level of the sampling hierarchy.}
  \item{relative}{Logical, if \code{TRUE} then alpha and beta diversity values are given relative to the value of gamma. This can be useful when comparing different indices.}
  \item{nsimul}{Number of permutation to use if \code{matr} is not of class 'permat'.}
  \item{control}{A list of arguments passed to quantitative permutation
          algorithms. If missing, the function  'permat.control' is used.}
  \item{x}{An object of class 'adipart'.}
  \item{\dots}{Other arguments passed to functions, e.g. base of logarithm for Shannon diversity (see \code{\link{diversity}}).}
}
\details{
Additive diversity partitioning means that mean alpha and beta diversity adds up to gamma diversity, thus beta diversity is measured in the same dimensions as alpha and gamma (Lande 1996). This additive procedure is than extended across multiple scales in a hierarchical sampling design with \eqn{i = 1, 2, 3, \ldots, m} levels of sampling (Crist et al. 2003). Samples in lower hierarchical levels are nested within higher level units, thus from \eqn{i=1} to \eqn{i=m} grain size is increasing under constant survey extent. At each level \eqn{i}, \eqn{\alpha_i} denotes average diversity found within samples.

At the highest sampling level, the diversity components are calculated as 

\deqn{\beta_m = \gamma - \alpha_m}{beta_m = gamma - alpha_m} 

For each lower sampling level as
\deqn{\beta_i = \alpha_{i+1} - \alpha_i}{beta_i = alpha_i+1 - alpha_i}

Then, the additive partition of diversity is 

\deqn{\gamma = \alpha_1 + \sum_{i=1}^m \beta_i}{gamma = alpha_1 + sum(beta_i)}

Average alpha components can be weighted uniformly (\code{weight="unif"}) to calculate it as simple average, or proportionally to sample abundances (\code{weight="prop"}) to calculate it as weighted average as follows

\deqn{\alpha_i = \sum_{j=1}^{n_i} D_{ij} w_{ij}}{alpha_i = sum(D_ij*w_ij)}

where \eqn{D_{ij}} is the diversity index and \eqn{w_{ij}} is the weight calculated for the \eqn{j}th sample at the \eqn{i}th sampling level.

The implementation of additive diversity partitioning follows Crist et al. 2003. It is based on species richness (\eqn{S}, not \eqn{S-1}), Shannon's and Simpson's diversity indices.

The expected diversity components are calculated \code{nsimul} times by individual based randomisation of the community data matrix. This is done by the \code{\link{permatfull}} and \code{\link{permatswap}} functions, and properties of the null model can be set by the \code{control} argument (see \code{\link{permat.control}}). The null matrics then evaluated via the function \code{\link{oecosimu}}.
}
\value{
An object of class 'adipart' with same structure as 'oecosimu' objects
}
\references{
Lande, R. 1996. Statistics and partitioning of species diversity, and similarity among multiple communities. \emph{Oikos}, 76, 5-13.

Crist, T.O., Veech, J.A., Gering, J.C. and Summerville, K.S. 2003. Partitioning species diversity across landscapes and regions: a hierarchical analysis of $\alpha$, $\beta$, and $\gamma$-diversity. \emph{Am. Nat.}, 162, 734-743.
}
\author{\enc{P\'eter S\'olymos}{Peter Solymos}, \email{solymos@ualberta.ca}}
\seealso{See \code{\link{permatfull}}, \code{\link{permatswap}} and \code{\link{permat.control}} for permutation settings, and \code{\link{oecosimu}} for calculating confidence levels.}
\examples{
data(mite)
data(mite.xy)
data(mite.env)
## Function to get equal area partitions of the mite data
cutter <- function (x, cut = seq(0, 10, by = 2.5)) {
    out <- rep(1, length(x))
    for (i in 2:(length(cut) - 1))
        out[which(x > cut[i] & x <= cut[(i + 1)])] <- i
    return(out)}
## The hierarchy of sample aggregation
levsm <- data.frame(
    l1=1:nrow(mite),
    l2=cutter(mite.xy$y, cut = seq(0, 10, by = 2.5)),
    l3=cutter(mite.xy$y, cut = seq(0, 10, by = 5)),
    l4=cutter(mite.xy$y, cut = seq(0, 10, by = 10)))
## Let's see in a map
par(mfrow=c(1,3))
plot(mite.xy, main="l1", col=levsm$l1+1)
plot(mite.xy, main="l2", col=levsm$l2+1)
plot(mite.xy, main="l3", col=levsm$l3+1)
par(mfrow=c(1,1))
## Additive diversity partitioning
adpMite <- adipart(mite, levsm, index="richness", nsimul=20)
adpMite
## Simple artificial example
set.seed(4321)
matr <- r2dtable(1, c(3,4,3,7,4,8,7,6), 3:9)[[1]]
strata <- cbind(1:8, rep(1:4,each=2), rep(1,8))
## restricted permutation within habitat classes
contr <- permat.control(hab = c("a","a","a","b","a","b","b","b"))
## Additive diversity partitioning ('oneway')
x1 <- adipart(matr, strata, index="shannon", nsimul=25, control=contr)
x1
}
\keyword{multivariate}

\encoding{UTF-8}
\name{adipart}
\alias{adipart}
\alias{plot.adipart}
\alias{print.adipart}
\alias{summary.adipart}
\alias{print.summary.adipart}
\title{Additive Diversity Partitioning}
\description{
In additive diversity partitioning, mean values of alpha diversity at lower levels of a sampling hierarchy are compared to the total number of species in the entire data set (gamma diversity) in the form: gamma = mean(alpha) + beta. Thus beta = gamma - mean(alpha).
}
\usage{
adipart(matr, strata, hclass = NULL, method="trad", 
index=c("richness", "shannon", "simpson"),
scales=seq(0, 2, 0.2), weights = "unif", test = TRUE, 
permtype = "full", times = 100, crit = 0.05, 
burnin = 10000, results = FALSE, ...)
\method{print}{adipart}(x, ...)
\method{summary}{adipart}(object, digits = 3, ...)
\method{print}{summary.adipart}(x, ...)
\method{plot}{adipart}(x, rel.yax = NULL, ymax = NULL, p.legend = "bottomright", ...)
}
\arguments{
  \item{matr}{a community data matrix with samples as rows and species as column, or a permutation object of class 'permat' (see Details).}
  \item{strata}{a matrix or data frame containing levels of sampling hierarchy as columns, columns from left to right will be treated as nested.}
  \item{hclass}{vector for coding habitat classes.}
  \item{method}{character, either \code{"trad"} for traditional diversity indices (see \code{index} argument for specifications), or \code{"tsallis"} for Tsallis diversity (see \code{scales} argument for specifications).}
  \item{index}{character, one of (if \code{habitat} is not \code{NULL}) or combination of (if \code{habitat = NULL}) \code{c("richness", "shannon", "simpson")}.}
  \item{scales}{vector for scales of the Tsallis diversity (see \code{\link{tsallis}} for details).}
  \item{weights}{character, \code{"unif"} for uniform weights, \code{"prop"} for weighting proportional to sample abundances to use in weighted averaging of individual alpha values within strata of a given level of the sampling hierarchy.}
  \item{test}{logical, whether a permutation test should be applied. If \ code{FALSE}, only observed values are returned.}
  \item{permtype}{character, \code{"full"} for permutation of community matrix via \code{\link{permatfull}}, \code{"swap"} for permutation via \code{\link{permatswap}}. Only used if \code{test = TRUE} and \code{matr} is not of class 'permat'.}
  \item{times}{number of permutation to use if \code{matr} is not of class 'permat'.}
  \item{crit}{two sided critical level for calculating confidence limits of the expected values. Only used if \code{test = TRUE}.}
  \item{burnin}{number of burnin steps when \code{permtype = "swap"} (see \code{\link{permatswap}}).}
  \item{results}{logical, whether null model results (individual alpha and beta diversity values) should be returned, only available for "oneway" design (\code{habitat = NULL}) if \code{test = TRUE}.}
  \item{x, object}{an object of class 'adipart'.}
  \item{digits}{number of significant digits to use in the output.}
  \item{rel.yax}{logical or \code{NULL}, \code{TRUE} for relative scaling of the y axis. The default (\code{NULL}) sets the scaling according to specific design ("oneway" or "twoway") criteria.}
  \item{ymax}{maximum of the vertical axis of the plot. Useful when \code{rel.yax = FALSE} and more than one plot should have same scales. \code{ymax} is ignored if it is lower than highest value to be plotted.}
  \item{p.legend}{position of the legend box for p-values, can be character of two coordinates, see \code{\link{legend}} for details. Use extreme large values to move the box out of the plotting region.}
  \item{\dots}{other arguments, e.g. arguments for \code{\link{permatfull}}, \code{\link{permatswap}}, \code{\link{par}} or \code{\link{print}}.}
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

The implementation of 'traditional' (\code{method="trad"}) 'oneway' (\code{hclass=NULL}) additive diversity partitioning follows Crist et al. 2003. This is 'traditional' in the sense, that it it is based on species richness (\eqn{S}, not \eqn{S-1}), Shannon's and Simpson's diversity indices. This is 'oneway' in the sense that habitat differences are assumed negligible (or, at least, well balanced or randomly sampled within strata), thus diversity partitioning is made according to successively larger grain sizes in the hierarchically nested sampling design.

The expected diversity components are calculated \code{times} times by individual based randomisation of the community data matrix. This is done by the \code{\link{permatfull}} (\code{permtype="full"}) or the \code{\link{permatswap}} (\code{permtype="swap"}) functions. Row and column sums are fixed by default, to change this behaviour, use the \code{fixedmar} argument. Restricted permutations can be set via the \code{reg} and \code{hab} arguments of these functions. Input objects can be either matrices, or objects of class 'permat'. For large community data sets and several thousands of permutations, it is advisable to use a matrix as input object to overcome memory usage problems (randomisation is made internally). If identical random matrices are needed for different computations, than it can be useful to make an object of 'permat' prior these analyses and use the 'permat' object as input.

The 'twoway' design (\code{hclass!=NULL}) can be used to compare diversity partitions among discrete habitat classes, and to calculate differentiation (beta diversity) among these habitat classes (Wagner et al. 2000). Within each habitat classes, the same sampling hierarchy is used, and among haitat diverity is calculated as the highest level beta component (\eqn{\beta_m}). The current implementation of null model testing for this 'twoway' design is based on the comparison of observed diversity components in a given habitat class with the expected components for all habitat classes.

The 'non-traditional' way of additive diversity partitioning is made via the Tsallis generalised entropy function (\code{method="tsallis"})  Scales of the scale parameter \eqn{q} can be set by the \code{scales} argument.
}
\note{
Please ensure that \code{strata} and \code{hclass} are meaningfully compiled (i.e. there are no missing states of the combinations), because at the moment there are no checks for that in the code. If number of observations within strata are very few, or not well balanced, than the permutation algorithm may fail to do acceptable randomisation.
}
\value{
An object of class 'adipart' with components:
  \item{input}{input objects (\code{m = matr}, \code{f = strata}, \code{h = habitat}) and the function call (\code{call}).}
  \item{obs}{observed diversity components (mean alpha, beta, and standard error for alpha).}
  \item{exp}{expected diversity components, both elements with items \code{p.value}, \code{mean}, lower and upper confidence limits (\code{cl1}, \code{cl2}) and starndardized effect size (\code{ses}, (observed - mean(expected)) / sd(expexted)):

    \code{alpha} expected alpha components,

    \code{beta} expected beta components.}
  \item{res}{null model distribution for the tested alpha and beta diversity values (\code{NULL} if \code{results = FALSE}):

    \code{alpha} null model distribution for alpha diversity, rows are permutations, columns are elements according to the \code{obs$alpha} matrix without the last row for gamma diversity,

    \code{beta} null model distribution for beta diversity, rows are permutations, columns are elements according to the \code{obs$beta} matrix.}
}
\references{
Lande, R. 1996. Statistics and partitioning of species diversity, and similarity among multiple communities. \emph{Oikos}, 76, 5-13.

Crist, T.O., Veech, J.A., Gering, J.C. and Summerville, K.S. 2003. Partitioning species diversity across landscapes and regions: a hierarchical analysis of $\alpha$, $\beta$, and $\gamma$-diversity. \emph{Am. Nat.}, 162, 734-743.

Wagner, H. H., Wildi, O. and Ewald, K.C. 2000. Additive partitioning of plant species diversity in an agricultural mosaic landscape. \emph{Landscape Ecology}, 15, 219-227.
}
\author{\enc{P\'eter S\'olymos}{Peter Solymos}, \email{solymos@ualberta.ca}}
\seealso{See \code{\link{permatfull}} and \code{\link{permatswap}} for permutation settings, and \code{\link{tsallis}} for Tsallis entropy.}
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
    l1=cutter(mite.xy$y, cut = seq(0, 10, by = 10)),
    l2=cutter(mite.xy$y, cut = seq(0, 10, by = 5)),
    l3=cutter(mite.xy$y, cut = seq(0, 10, by = 2.5)),
    l4=1:nrow(mite))
## Let's see in a map
par(mfrow=c(1,3))
plot(mite.xy, main="l1", col=levsm$l1+1)
plot(mite.xy, main="l2", col=levsm$l2+1)
plot(mite.xy, main="l3", col=levsm$l3+1)
par(mfrow=c(1,1))
## Additive diversity partitioning
adpMite <- adipart(mite, levsm, index=c("ri", "sh", "si"), times=10)
adpMite
summary(adpMite)
plot(adpMite)
## Simple artificial example
set.seed(4321)
matr <- r2dtable(1, c(3,4,3,7,4,8,7,6), 3:9)[[1]]
habitat <- c("a","a","a","b","a","b","b","b")
strata <- cbind(
   rep(1,8),
   rep(1:4,each=2),
   1:8)
## Additive diversity partitioning ('oneway')
x1 <- adipart(matr, strata, index=c("ri", "sh", "si"), times=25)
x1
summary(x1)
plot(x1)
## 'Twoway' partitioning
x2 <- adipart(matr, strata, habitat, times=25)
x2
summary(x2)
plot(x2)
}
\keyword{multivariate}

\encoding{UTF-8}
\name{adipart}
\alias{adipart}
\alias{multipart}
\alias{print.adipart}
\alias{print.multipart}
\alias{hiersimu}
\alias{print.hiersimu}
\title{Additive Diversity Partitioning and Hierarchical Null Model Testing}
\description{
In diversity partitioning, mean values of alpha diversity at lower levels of a sampling 
hierarchy are compared to the total diversity in the entire data set (gamma diversity). 
The prtitioning can be additive (\code{adipart}) or multiplicative (\code{multipart}). 
In hierarchical null model testing, a statistic returned by a function is evaluated 
according to a nested hierarchical sampling design (\code{hiersimu}).
}
\usage{
adipart(formula, data, index=c("richness", "shannon", "simpson"),
    weights=c("unif", "prop"), relative = FALSE, nsimul=99, control, ...)
multipart(formula, data, index=c("renyi", "tsallis"), scales = 1,
    global = FALSE, relative = FALSE, nsimul=99, control, ...)
hiersimu(formula, data, FUN, location = c("mean", "median"),
    relative = FALSE, drop.highest = FALSE, nsimul=99, control, ...)
\method{print}{adipart}(x, ...)
\method{print}{multipart}(x, ...)
\method{print}{hiersimu}(x, ...)
}
\arguments{
  \item{formula}{A two sided model formula in the form \code{y ~ x}, where \code{y} 
    is the community data matrix with samples as rows and species as column. Right 
    hand side (\code{x}) must contain factors referring to levels of sampling hierarchy, 
    terms from right to left will be treated as nested (first column is the lowest, 
    last is the highest level). These variables must be factors in order to unambiguous 
    handling. Interaction terms are not allowed.}
  \item{data}{A data frame where to look for variables defined in the right hand side 
    of \code{formula}. If missing, variables are looked in the global environment.}
  \item{index}{Character, the diversity index to be calculated (see Details).}
  \item{weights}{Character, \code{"unif"} for uniform weights, \code{"prop"} for 
    weighting proportional to sample abundances to use in weighted averaging of individual 
    alpha values within strata of a given level of the sampling hierarchy.}
  \item{relative}{Logical, if \code{TRUE} then alpha and beta diversity values are given 
    relative to the value of gamma for function \code{adipart}. For function \code{multipart},
    it sets the standardization of the beta diversity value with its maximum (see Details).}
  \item{scales}{Numeroc, of length 1, the order of the generalized diversity index 
    to be used.}
  \item{global}{Logical, indicates the calculation of beta values, see Details.}
  \item{nsimul}{Number of permutation to use if \code{matr} is not of class 'permat'.
    If \code{nsimul = 0}, only the \code{FUN} argument is evaluated. It is thus possible
    to reuse the statistic values without using a null model.}
  \item{control}{A list of arguments passed to quantitative permutation
          algorithms. If missing, the function  'permat.control' is used.}
  \item{FUN}{A function to be used by \code{hiersimu}.}
  \item{location}{Character, identifies which function (mean or median) is used to 
    calculate location of the samples.}
  \item{drop.highest}{Logical, to drop the highest level or not. When \code{FUN} 
    evaluates only arrays with at least 2 dimensions, highest level should be dropped, 
    or not selected at all.}
  \item{x}{An object to print.}
  \item{\dots}{Other arguments passed to functions, e.g. base of logarithm for 
    Shannon diversity (see \code{\link{diversity}}).}
}
\details{
Additive diversity partitioning means that mean alpha and beta diversity adds up to gamma 
diversity, thus beta diversity is measured in the same dimensions as alpha and gamma 
(Lande 1996). This additive procedure is than extended across multiple scales in a 
hierarchical sampling design with \eqn{i = 1, 2, 3, \ldots, m} levels of sampling 
(Crist et al. 2003). Samples in lower hierarchical levels are nested within higher level 
units, thus from \eqn{i=1} to \eqn{i=m} grain size is increasing under constant survey 
extent. At each level \eqn{i}, \eqn{\alpha_i} denotes average diversity found within samples.

At the highest sampling level, the diversity components are calculated as 
\deqn{\beta_m = \gamma - \alpha_m}{beta_m = gamma - alpha_m} 
For each lower sampling level as
\deqn{\beta_i = \alpha_{i+1} - \alpha_i}{beta_i = alpha_i+1 - alpha_i}
Then, the additive partition of diversity is 
\deqn{\gamma = \alpha_1 + \sum_{i=1}^m \beta_i}{gamma = alpha_1 + sum(beta_i)}

Average alpha components can be weighted uniformly (\code{weight="unif"}) to calculate 
it as simple average, or proportionally to sample abundances (\code{weight="prop"}) to 
calculate it as weighted average as follows
\deqn{\alpha_i = \sum_{j=1}^{n_i} D_{ij} w_{ij}}{alpha_i = sum(D_ij*w_ij)}
where \eqn{D_{ij}} is the diversity index and \eqn{w_{ij}} is the weight calculated for 
the \eqn{j}th sample at the \eqn{i}th sampling level.

The implementation of additive diversity partitioning in \code{adipart} follows Crist et 
al. 2003. It is based on species richness (\eqn{S}, not \eqn{S-1}), Shannon's and 
Simpson's diversity indices states as argument \code{index}.

Multiplicative diversity partitioning is based on Whittaker's (1972) ideas, that has 
recently been generalised to one parametric diversity families (i.e. \enc{R\'enyi}{Renyi} 
and Tsallis) by Jost (2006, 2007). Jost recomments to use the numbers equivalents 
(Hill numbers), instead of pure diversities, and proofs, that this satisfies the 
multiplicative partitioning requirements.

The current implementation of \code{multipart} calculates Hill numbers based on the 
functions \code{\link{renyi}} and \code{\link{tsallis}} (provided as \code{index} argument). 
If values for more than one \code{scales} values are desired, it should be done in separate 
runs, because it adds extra dimensionality to the implementation, which has not been resolved 
efficiently.

Alpha diversities are then the averages of these Hill numbers for each hierarchy levels, 
the global gamma diversity is the alpha value calculated for the highest hierarchy level. 
When \code{global = TRUE}, beta is calculated relative to the global gamma value:
\deqn{\beta_i = \gamma / \alpha_{i}}{beta_i = gamma / alpha_i}
when \code{global = FALSE], beta is calculated relative to local gamma values (local gamma
means the diversity the calculated for a particular cluster based on the pooled abundance vector):
\deqn{\beta_ij = \alpha_{(i+1)j} / mean(\alpha_{ij})}{beta_ij = alpha_(i+1)j / mean(alpha_i)}
where \eqn{j} is a particular cluster at hierarchy level \eqn{i}. Then beta diversity value for
level \eqn{i} is the mean of the beta values of the clusters at that level,
\eqn{\beta_{i} = mean(\beta_{ij})}.

If \code{relative = TRUE} for \code{multipart}, the respective beta diversity values are
standardized by their maximum expected values (\eqn{mean(\beta_{ij}) / \beta_{max,ij}})
given as \eqn{\beta_{max,ij} = n_{j}} (the number of lower level units in a given cluster \eqn{j}).

The expected diversity components are calculated \code{nsimul} times by individual based 
randomisation of the community data matrix. This is done by the \code{\link{permatfull}} 
and \code{\link{permatswap}} functions, and properties of the null model can be set by 
the \code{control} argument (see \code{\link{permat.control}}). The null matrics then 
evaluated via the function \code{\link{oecosimu}}.

\code{hiersimu} works almost the same as \code{adipart}, but without comparing the actual 
statistic values returned by \code{FUN} to the highest possible value (cf. gamma diversity). 
This is so, because in most of the cases, it is difficult to ensure additive properties of 
the mean statistic values along the hierarchy.
}
\value{
An object of class 'adipart' or 'hiersimu' with same structure as 'oecosimu' objects.
}

\references{
Crist, T.O., Veech, J.A., Gering, J.C. and Summerville,
K.S. (2003). Partitioning species diversity across landscapes and regions:
a hierarchical analysis of \eqn{\alpha}, \eqn{\beta}, and
\eqn{\gamma}-diversity.
\emph{Am. Nat.}, \bold{162}, 734--743.

Jost, L. (2006). Entropy and diversity.
\emph{Oikos}, \bold{113}, 363--375.

Jost, L. (2007). Partitioning diversity into independent alpha and beta components.
\emph{Ecology}, \bold{88}, 2427--2439.

Lande, R. (1996). Statistics and partitioning of species
diversity, and similarity among multiple communities.
\emph{Oikos}, \bold{76}, 5--13.

Whittaker, R. (1972). Evolution and measurement of species diversity.
\emph{Taxon}, \bold{21}, 213--251.
}

\author{\enc{P\'eter S\'olymos}{Peter Solymos}, \email{solymos@ualberta.ca}}
\seealso{See \code{\link{permatfull}}, \code{\link{permatswap}} and \code{\link{permat.control}} 
  for permutation settings, and \code{\link{oecosimu}} for calculating confidence levels.}
\examples{
data(mite)
data(mite.xy)
data(mite.env)
## Function to get equal area partitions of the mite data
cutter <- function (x, cut = seq(0, 10, by = 2.5)) {
    out <- rep(1, length(x))
    for (i in 2:(length(cut) - 1))
        out[which(x > cut[i] & x <= cut[(i + 1)])] <- i
    return(as.factor(out))}
## The hierarchy of sample aggregation
levsm <- data.frame(
    l1=as.factor(1:nrow(mite)),
    l2=cutter(mite.xy$y, cut = seq(0, 10, by = 2.5)),
    l3=cutter(mite.xy$y, cut = seq(0, 10, by = 5)),
    l4=cutter(mite.xy$y, cut = seq(0, 10, by = 10)))
## Let's see in a map
par(mfrow=c(1,3))
plot(mite.xy, main="l1", col=as.numeric(levsm$l1)+1)
plot(mite.xy, main="l2", col=as.numeric(levsm$l2)+1)
plot(mite.xy, main="l3", col=as.numeric(levsm$l3)+1)
par(mfrow=c(1,1))
## Additive diversity partitioning
adpMite <- adipart(mite ~., levsm, index="richness", nsimul=20)
adpMite
## Simple artificial example
set.seed(4321)
matr <- r2dtable(1, c(3,4,3,7,4,8,7,6), 3:9)[[1]]
strata <- data.frame(letters[1:8], rep(letters[1:4],each=2), rep("a",8))
## Restricted permutation within habitat classes
contr <- permat.control(strata = c("a","a","a","b","a","b","b","b"))
## Additive diversity partitioning
x1 <- adipart(matr ~ ., strata, index="shannon", nsimul=25, control=contr)
x1
## Multiplicative diversity partitioning
x2 <- multipart(matr ~ ., strata, index="renyi", scales=1, nsimul=25)
x2
## Hierarchical null model testing
## diversity analysis (similar to adipart)
hiersimu(mite ~., levsm, diversity, relative=TRUE, nsimul=25)
## Hierarchical testing with the Morisita index
morfun <- function(x) dispindmorisita(x)$imst
hiersimu(mite ~., levsm, morfun, drop.highest=TRUE, nsimul=25)
}
\keyword{multivariate}

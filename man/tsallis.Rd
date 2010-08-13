\encoding{UTF-8}
\name{tsallis}
\alias{tsallis}
\alias{tsallisaccum}
\alias{persp.tsallisaccum}
\title{Tsallis Diversity and Corresponding Accumulation Curves}
\description{
Function \code{tsallis} find Tsallis diversities with any scale or the corresponding evenness measures. Function \code{tsallisaccum} finds these statistics with accumulating sites.
}
\usage{
tsallis(x, scales = seq(0, 2, 0.2), norm = FALSE, hill = FALSE)
tsallisaccum(x, scales = seq(0, 2, 0.2), permutations = 100, raw = FALSE, ...)
\method{persp}{tsallisaccum}(x, theta = 220, phi = 15, col = heat.colors(100), zlim, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Community data matrix or plotting object. }
  \item{scales}{Scales of Tsallis diversity.}
  \item{norm}{Logical, if \code{TRUE} diversity values are normalized by their maximum (diversity value at equiprobability conditions).}
  \item{hill}{Calculate Hill numbers.}
  \item{permutations}{Number of random permutations in accumulating sites.}
  \item{raw}{If \code{FALSE} then return summary statistics of permutations, and if TRUE then returns the individual permutations.}
  \item{theta, phi}{angles defining the viewing direction. \code{theta} gives the azimuthal direction and \code{phi} the colatitude.}
  \item{col}{Colours used for surface.}
  \item{zlim}{Limits of vertical axis.}
  \item{\dots}{Other arguments which are passed to \code{tsallis} and to graphical functions.}

}
\details{
The Tsallis diversity (also equivalent to Patil and Taillie diversity) is a one-parametric generalised entropy function, defined as:

\deqn{H_q = \frac{1}{q-1} (1-\sum_{i=1}^S p_i^q)}{H.q = 1/(q-1)(1-sum(p^q))}

where \eqn{q} is a scale parameter, \eqn{S} the number of species in the sample (Tsallis 1988, Tothmeresz 1995). This diversity is concave for all \eqn{q>0}, but non-additive (Keylock 2005). For \eqn{q=0} it gives the number of species minus one, as \eqn{q} tends to 1 this gives Shannon diversity, for \eqn{q=2} this gives the Simpson index (see function \code{\link{diversity}}).

If \code{norm = TRUE}, \code{tsallis} gives values normalized by the maximum:

\deqn{H_q(max) = \frac{S^{1-q}-1}{1-q}}{H.q(max) = (S^(1-q)-1)/(1-q)}

where \eqn{S} is the number of species. As \eqn{q} tends to 1, maximum is defined as \eqn{ln(S)}.

If \code{hill = TRUE}, \code{tsallis} gives Hill numbers (numbers equivalents, see Jost 2007):

\deqn{D_q = (1-(q-1) H)^{1/(1-q)}}{D.q = (1-(q-1)*H)^(1/(1-q))}

Details on plotting methods and accumulating values can be found on the help pages of the functions \code{\link{renyi}} and \code{\link{renyiaccum}}.
}
\value{
Function \code{tsallis} returns a data frame of selected indices. Function \code{tsallisaccum} with argument \code{raw = FALSE} returns a three-dimensional array, where the first dimension are the accumulated sites, second dimension are the diversity scales, and third dimension are the summary statistics \code{mean}, \code{stdev}, \code{min}, \code{max}, \code{Qnt 0.025} and \code{Qnt 0.975}. With argument \code{raw = TRUE} the statistics on the third dimension are replaced with individual permutation results.
}
\references{
Tsallis, C. (1988) Possible generalization of Boltzmann-Gibbs statistics. 
  \emph{J. Stat. Phis.} 52, 479--487.

Tothmeresz, B. (1995) Comparison of different methods for diversity
  ordering. \emph{Journal of Vegetation Science} \bold{6}, 283--290.

Patil, G. P. and Taillie, C. (1982) Diversity as a concept and its measurement.
  \emph{J. Am. Stat. Ass.} \bold{77}, 548--567.

Keylock, C. J. (2005) Simpson diversity and the Shannon-Wiener index as special cases of a generalized entropy.
  \emph{Oikos} \bold{109}, 203--207.

Jost, L (2007) Partitioning diversity into independent alpha and beta components.
  \emph{Ecology} \bold{88}, 2427--2439.
}
\author{\enc{Péter Sólymos}{Peter Solymos}, \email{solymos@ualberta.ca}, based on the code of Roeland Kindt and Jari Oksanen written for \code{renyi}}
\seealso{
Plotting methods and accumulation routines are based on functions \code{\link{renyi}} and \code{\link{renyiaccum}}. An object of class 'tsallisaccum' can be used with function \code{\link{rgl.renyiaccum}} as well. See also settings for \code{\link{persp}}.
}
\examples{
data(BCI)
i <- sample(nrow(BCI), 12)
x1 <- tsallis(BCI[i,])
x1
diversity(BCI[i,],"simpson") == x1[["2"]]
plot(x1)
x2 <- tsallis(BCI[i,],norm=TRUE)
x2
plot(x2)
mod1 <- tsallisaccum(BCI[i,])
plot(mod1, as.table=TRUE, col = c(1, 2, 2))
persp(mod1)
mod2 <- tsallisaccum(BCI[i,], norm=TRUE)
persp(mod2,theta=100,phi=30)
}
\keyword{multivariate}

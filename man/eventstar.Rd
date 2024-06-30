\encoding{UTF-8}
\name{eventstar}
\alias{eventstar}

\title{
Scale Parameter at the Minimum of the Tsallis Evenness Profile
}
\description{
The function \code{eventstar} finds the minimum (\eqn{q^*}{q*}) of the 
evenness profile based on the Tsallis entropy. This scale factor
of the entropy represents a specific weighting of species
relative frequencies that leads to minimum evenness of the
community (Mendes et al. 2008).
}
\usage{
eventstar(x, qmax = 5)
}
\arguments{
  \item{x}{
A community matrix or a numeric vector.
}
  \item{qmax}{
Maximum scale parameter of the Tsallis entropy to be used in 
finding the minimum of Tsallis based evenness
in the range \code{c(0, qmax)}.
}
}
\details{
The function \code{eventstar} finds a characteristic value of the scale 
parameter \eqn{q} of the Tsallis entropy corresponding to
minimum of the evenness (equitability) profile based on Tsallis entropy.
This value was proposed by Mendes et al. (2008) as \eqn{q^*}{q*}.

The \eqn{q^\ast}{q*} index represents the scale parameter of
the one parameter Tsallis diversity family that leads to
the greatest deviation from the maximum equitability given the relative 
abundance vector of a community.

The value of \eqn{q^\ast}{q*} is found by identifying the minimum
of the evenness profile over scaling factor \eqn{q}{q} by
one-dimensional minimization. Because evenness profile is
known to be a convex function, it is guaranteed that underlying
\code{\link{optimize}} function will find a unique solution
if it is in the range \code{c(0, qmax)}.

The scale parameter value \eqn{q^\ast}{q*} is used to 
find corresponding values of diversity (\eqn{H_{q^\ast}}{H.q*}), 
evenness (\eqn{H_{q^\ast}(\max)}{H.q*(max)}),
and numbers equivalent (\eqn{D_{q^\ast}}{D.q*}). For calculation
details, see \code{\link{tsallis}} and Examples below.

Mendes et al. (2008) advocated the use of \eqn{q^\ast}{q*}
and corresponding diversity, evenness, and Hill numbers, because
it is a unique value representing the diversity profile, and is
is positively associated with rare species in the community,
thus it is a potentially useful indicator of certain
relative abundance distributions of the communities.
}
\value{
A data frame with columns:

  \item{\code{qstar}}{ scale parameter value \eqn{q\ast}{q*}
    corresponding to minimum value of Tsallis based evenness profile.}
  \item{\code{Estar}}{ Value of evenness based on normalized Tsallis 
    entropy at \eqn{q^\ast}{q*}.}
  \item{\code{Hstar}}{ Value of Tsallis entropy at \eqn{q^\ast}{q*}.}
  \item{\code{Dstar}}{ Value of Tsallis entropy at \eqn{q^\ast}{q*} 
    converted to numbers equivalents
    (also called as Hill numbers, effective number of species, 
    \sQuote{true} diversity; cf. Jost 2007).}

See \code{\link{tsallis}} for calculation details.
}
\references{
Mendes, R.S., Evangelista, L.R., Thomaz, S.M.,
  Agostinho, A.A. and Gomes, L.C. (2008) A unified
  index to measure ecological diversity and species
  rarity. \emph{Ecography} \bold{31}, 450--456.

Jost, L. (2007) Partitioning diversity into independent alpha and beta components.
  \emph{Ecology} \bold{88}, 2427--2439.

Tsallis, C. (1988) Possible generalization of Boltzmann-Gibbs statistics. 
  \emph{J. Stat. Phis.} \bold{52}, 479--487.
}
\note{
Values for \eqn{q^\ast}{q*} found by Mendes et al. (2008) ranged
from 0.56 and 1.12 presenting low variability, so an
interval between 0 and 5 should safely encompass
the possibly expected \eqn{q^\ast}{q*} values in practice,
but profiling the evenness and changing the value of
the \code{qmax} argument is advised if output values
near the range limits are found.
}
\author{
Eduardo Ribeiro Cunha \email{edurcunha@gmail.com} and 
Heloisa Beatriz Antoniazi Evangelista \email{helobeatriz@gmail.com}, 
with technical input of Péter Sólymos.
}
\seealso{
Tsallis entropy: \code{\link{tsallis}}
}
\examples{
data(BCI)
(x <- eventstar(BCI[1:5,]))
## profiling
y <- as.numeric(BCI[10,])
(z <- eventstar(y))
q <- seq(0, 2, 0.05)
Eprof <- tsallis(y, scales=q, norm=TRUE)
Hprof <- tsallis(y, scales=q)
Dprof <- tsallis(y, scales=q, hill=TRUE)
opar <- par(mfrow=c(3,1))
plot(q, Eprof, type="l", main="Evenness")
abline(v=z$qstar, h=tsallis(y, scales=z$qstar, norm=TRUE), col=2)
plot(q, Hprof, type="l", main="Diversity")
abline(v=z$qstar, h=tsallis(y, scales=z$qstar), col=2)
plot(q, Dprof, type="l", main="Effective number of species")
abline(v=z$qstar, h=tsallis(y, scales=z$qstar, hill=TRUE), col=2)
par(opar)
}
\keyword{ optimize }
\keyword{ multivariate }
\keyword{ utilities }

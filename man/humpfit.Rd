\name{humpfit}
\alias{humpfit}
\alias{summary.humpfit}
\alias{lines.humpfit}
\alias{plot.humpfit}
\alias{points.humpfit}
\alias{predict.humpfit}
\alias{profile.humpfit}

\title{No-interaction Model for Hump-backed Species Richness vs. Biomass }
\description{
  Function \code{humpfit} fits a no-interaction model for species
  richness vs. biomass data (Oksanen 1996). This is a null model that
  produces a hump-backed response as an artifact of plant size and
  density. 
}
\usage{
humpfit(mass, spno, family = poisson, start)
\method{summary}{humpfit}(object, ...)
\method{predict}{humpfit}(object, newdata = NULL, ...) 
\method{plot}{humpfit}(x, xlab = "Biomass", ylab = "Species Richness", lwd = 2, 
    l.col = "blue", p.col = 1, type = "b", ...)
\method{points}{humpfit}(x, ...)
\method{lines}{humpfit}(x, segments=101,  ...)
\method{profile}{humpfit}(fitted, parm = 1:3, alpha = 0.01, maxsteps = 20, del = zmax/5, ...)
}

\arguments{
  \item{mass}{Biomass. }
  \item{spno}{Species richness.}
  \item{start}{Vector of starting values for all three parameters.}
  \item{family}{Family of error distribution. Any \code{\link{family}}
    can be used, but the link function is always Fisher's diversity
    model, and other \code{link} functions are silently ignored. }
  \item{x, object, fitted}{Result object of \code{humpfit}}
  \item{newdata}{Values of \code{mass} used in \code{predict}. The
    original data values are used if missing.}
  \item{xlab,ylab}{Axis labels in \code{plot}}
  \item{lwd}{Line width}
  \item{l.col, p.col}{Line and point colour in \code{plot}}
  \item{type}{Type of \code{plot}: \code{"p"} for observed points,
    \code{"l"} for fitted lines, \code{"b"} for both, and \code{"n"} for
    only setting axes.}
  \item{segments}{Number of segments used for fitted lines.}
  \item{parm}{Profiled parameters.}
  \item{alpha, maxsteps, del}{Parameters for profiling range and density.}
  \item{...}{Other parameters to functions.}
}
\details{
  The no-interaction model assumes that the humped species richness
  pattern along biomass gradient is an artifact of plant size and
  density (Oksanen 1996). For low-biomass sites, it assumes that plants
  have a fixed size, and biomass increases with increasing number of
  plants. When the sites becomes crowded, the number of plants and
  species richness reaches the maximum. Higher biomass is reached by
  increasing the plant size, and then the number of plants and species
  richness will decrease. At biomasses below the hump, plant number and
  biomass are linearly related, and above the hump, plant number is
  proportional to inverse squared biomass. The number of plants is
  related to the number of species by the relationship (\code{link}
  function) from Fisher's log-series (Fisher et al. 1943).

  The parameters of the model are:
  \enumerate{
    \item \code{hump}: the location of the hump on the biomass gradient.
    \item \code{scale}: an arbitrary multiplier to translate the biomass
    into virtual number of plants.
    \item \code{alpha}: Fisher's \eqn{\alpha}{alpha} to translate the
    virtual number of plants into number of species.
  }
  The parameters \code{scale} and \code{alpha} are intermingled and this
  function should not be used for estimating Fisher's
  \eqn{\alpha}{alpha}.  Probably the only meaningful and interesting
  parameter is the location of the \code{hump}.

  Function may be very
  difficult to fit and easily gets trapped into local solutions, or
  fails with non-Poisson families, and function \code{profile} should be
  used to inspect the fitted models. If you have loaded \code{package}
  \pkg{MASS}, you can use functions
  \code{\link[MASS]{plot.profile}},
  \code{\link[MASS]{pairs.profile}} for graphical inspection of the
  profiles, and \code{\link[MASS]{confint.profile.glm}} for the profile
  based confidence intervals. 

  The original model intended to show that there is no need to
  speculate about \dQuote{competition} and \dQuote{stress}
  (Al-Mufti et al. 1977),
  but humped response can be produced as an artifact of using fixed
  plot size for varying plant sizes and densities.  } \value{ The
  function returns an object of class \code{"humpfit"} inheriting from
  class \code{"glm"}. The result object has specific \code{summary},
  \code{predict}, \code{plot}, \code{points} and \code{lines}
  methods. In addition, it can be accessed by the following methods
  for \code{\link{glm}} objects: \code{\link{AIC}},
  \code{\link{extractAIC}}, \code{\link{deviance}},
  \code{\link{coef}}, \code{\link{residuals.glm}} (except 
  \code{type = "partial"}), \code{\link{fitted}}, and perhaps some others. 
  In addition, function \code{\link[ellipse]{ellipse.glm}} (package
  \pkg{ellipse}) can be used to draw approximate confidence ellipses
  for pairs of parameters, if the normal assumptions look appropriate.
  } 
\references{ 
  Al-Mufti, M.M., Sykes, C.L, Furness, S.B., Grime, J.P
  & Band, S.R. (1977) A quantitative analysis of shoot phenology and
  dominance in herbaceous vegetation. \emph{Journal of Ecology}
  65,759--791.
  
  Fisher, R.A., Corbet, A.S. & Williams, C.B. (1943) The relation
  between the number of species and the number of individuals in a
  random sample of of an animal population. \emph{Journal of Animal
    Ecology} 12, 42--58.
  
  Oksanen, J. (1996) Is the humped relationship between species richness
  and biomass an artefact due to plot size? \emph{Journal of Ecology}
  84, 293--295.
}
\author{ Jari Oksanen }
\note{
  The function is a replacement for the original \code{GLIM4} function
  at the archive of Journal of Ecology.  There the function was
  represented as a mixed \code{\link{glm}} with one non-linear
  parameter (\code{hump}) and a special one-parameter link function from
  Fisher's log-series.  The current function directly applies non-linear
  maximum likelihood fitting using function \code{\link{nlm}}.  Some
  expected problems with the current approach are:
  \itemize{
    \item The function is discontinuous at \code{hump} and may be
    difficult to optimize in some cases (the lines will always join, but
    the derivative jumps).
    \item The function does not try very hard to find sensible starting
    values and can fail. The user may supply starting values in
    argument \code{start} if fitting fails.
    \item The estimation is unconstrained, but both \code{scale} and
    \code{alpha} should always be positive.  Perhaps they should be
    fitted as logarithmic. Fitting \code{\link{Gamma}}
    \code{\link{family}} models might become easier, too.
    }
}


\seealso{\code{\link{fisherfit}}, \code{\link[MASS]{profile.glm}},
  \code{\link[MASS]{confint.glm}}. }
\examples{
##
## Data approximated from Al-Mufti et al. (1977)
##
\donttest{
mass <- c(140,230,310,310,400,510,610,670,860,900,1050,1160,1900,2480)
spno <- c(1,  4,  3,  9, 18, 30, 20, 14,  3,  2,  3,  2,  5,  2)
sol <- humpfit(mass, spno)
summary(sol) # Almost infinite alpha...
plot(sol)
# confint is in MASS, and impicitly calls profile.humpfit.
# Parameter 3 (alpha) is too extreme for profile and confint, and we
# must use only "hump" and "scale".
library(MASS)
plot(profile(sol, parm=1:2))
confint(sol, parm=c(1,2))
}
}
\keyword{models }
\keyword{regression }
\keyword{nonlinear}

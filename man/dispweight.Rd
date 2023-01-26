\encoding{UTF-8}
\name{dispweight}
\alias{dispweight}
\alias{gdispweight}
\alias{summary.dispweight}
\title{Dispersion-based weighting of species counts}

\description{Transform abundance data downweighting species that are 
  overdispersed to the Poisson error.}

\usage{
dispweight(comm, groups, nsimul = 999, nullmodel = "c0_ind",
    plimit = 0.05)
gdispweight(formula, data, plimit = 0.05)
\method{summary}{dispweight}(object, ...)
}

\arguments{
  \item{comm}{Community data matrix.}
  \item{groups}{Factor describing the group structure. If missing, all 
     sites are regarded as belonging to one group. \code{NA} values are 
     not allowed.}
  \item{nsimul}{Number of simulations.}
  \item{nullmodel}{The \code{\link{nullmodel}} used in
    \code{\link{commsim}} within \code{groups}. The default
    follows Clarke et al. (2006).}
  \item{plimit}{Downweight species if their \eqn{p}-value is at or
     below this limit.}
  \item{formula, data}{Formula where the left-hand side is the
    community data frame and right-hand side gives the explanatory
    variables. The explanatory variables are found in the data frame
    given in \code{data} or in the parent frame.}
  \item{object}{Result object from \code{dispweight} or
    \code{gdispweight}.}
  \item{\dots}{Other parameters passed to functions.}
}

\details{

The dispersion index (\eqn{D}) is calculated as ratio between variance
and expected value for each species.  If the species abundances follow
Poisson distribution, expected dispersion is \eqn{E(D) = 1}, and if
\eqn{D > 1}, the species is overdispersed. The inverse \eqn{1/D} can
be used to downweight species abundances.  Species are only
downweighted when overdispersion is judged to be statistically
significant (Clarke et al. 2006).

Function \code{dispweight} implements the original procedure of Clarke
et al. (2006). Only one factor can be used to group the sites and to
find the species means. The significance of overdispersion is assessed
freely distributing individuals of each species within factor
levels. This is achieved by using \code{\link{nullmodel}}
\code{"c0_ind"} (which accords to Clarke et al. 2006), but other
nullmodels can be used, though they may not be meaningful (see
\code{\link{commsim}} for alternatives). If a species is absent in
some factor level, the whole level is ignored in calculation of
overdispersion, and the number of degrees of freedom can vary among
species. The reduced number of degrees of freedom is used as a divisor
for overdispersion \eqn{D}, and such species have higher dispersion
and hence lower weights in transformation.

Function \code{gdispweight} is a generalized parametric version of
\code{dispweight}. The function is based on \code{\link{glm}} with
\code{\link{quasipoisson}} error \code{\link{family}}. Any
\code{\link{glm}} model can be used, including several factors or
continuous covariates. Function \code{gdispweight} uses the same test
statistic as \code{dispweight} (Pearson Chi-square), but it does not
ignore factor levels where species is absent, and the number of
degrees of freedom is equal for all species. Therefore transformation
weights can be higher than in \code{dispweight}. The
\code{gdispweight} function evaluates the significance of
overdispersion parametrically from Chi-square distribution
(\code{\link{pchisq}}).

Functions \code{dispweight} and \code{gdispweight} transform data, but
they add information on overdispersion and weights as attributes of
the result. The \code{summary} can be used to extract and print that
information.  
}

\value{
Function returns transformed data with the following new attributes:
    \item{D}{Dispersion statistic.}
    \item{df}{Degrees of freedom for each species.}
    \item{p}{\eqn{p}-value of the Dispersion statistic \eqn{D}.}
    \item{weights}{weights applied to community data.}
    \item{nsimul}{Number of simulations used to assess the \eqn{p}-value,
      or \code{NA} when simulations were not performed.}
    \item{nullmodel}{The name of \code{\link{commsim}} null model, or
      \code{NA} when simulations were not performed.}
}

\references{
Clarke, K. R., M. G. Chapman, P. J. Somerfield, and
H. R. Needham. 2006. Dispersion-based weighting of species counts in
assemblage analyses. \emph{Marine Ecology Progress Series}, 320,
11–27.
}

\author{
  Eduard Szöcs \email{eduardszoesc@gmail.com} wrote the original
  \code{dispweight}, Jari Oksanen significantly modified the code,
  provided support functions and developed \code{gdispweight}.
}


\examples{
data(mite, mite.env)
## dispweight and its summary
mite.dw <- with(mite.env, dispweight(mite, Shrub, nsimul = 99))
## IGNORE_RDIFF_BEGIN
summary(mite.dw)
## IGNORE_RDIFF_END
## generalized dispersion weighting
mite.dw <- gdispweight(mite ~ Shrub + WatrCont, data = mite.env)
rda(mite.dw ~ Shrub + WatrCont, data = mite.env)
}


\keyword{multivariate}
\keyword{manip}


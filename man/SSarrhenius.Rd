\name{SSarrhenius}
\alias{SSarrhenius}
\alias{SSlomolino}
\alias{SSgitay}
\alias{SSgleason}

\title{
  Self-Starting nls Species-Area Models
}

\description{
  These functions provide self-starting species-area models for
  non-linear regression (\code{\link{nls}}). They can also be used for
  fitting species accumulation models in
  \code{\link{fitspecaccum}}. These models (and many more) are reviewed
  by Dengler (2009).
}

\usage{
SSarrhenius(area, k, z)
SSgleason(area, k, slope)
SSgitay(area, k, slope)
SSlomolino(area, Asym, xmid, slope)
}

\arguments{
  \item{area}{
    Area or size of the sample: the independent variable.
  }
  \item{k, z, slope, Asym, xmid}{
    Estimated model parameters: see Details.
  }
}
\details{
  All these functions are assumed to be used for species richness
  (number of species) as the independent variable, and area or sample
  size as the independent variable. Basically, these define least
  squares models of untransformed data, and will differ from models
  for transformed species richness or models with non-Gaussian error.

  The Arrhenius model (\code{SSarrhenius}) is the expression
  \code{k*area^z}. This is the most classical model that can be found in
  any textbook of ecology (and also in Dengler 2009). Parameter \code{z}
  is the steepness of the species-area curve, and \code{k} is the
  expected number of species in a unit area.

  The Gleason model (\code{SSgleason}) is a linear expression 
  \code{k + slope*log(area)} (Dengler 200). This is a linear model,  
  and starting values give the final estimates; it is provided to 
  ease comparison with other models.

  The Gitay model (\code{SSgitay}) is a quadratic logarithmic expression
  \code{(k + slope*log(area))^2} (Gitay et al. 1991, Dengler
  2009). Parameter \code{slope} is the steepness of the species-area
  curve, and \code{k} is the square root of expected richness in a unit
  area. 

  The Lomolino model (\code{SSlomolino}) is
  \code{Asym/(1 + slope^log(xmid/area))} (Lomolino 2000, Dengler 2009).
  Parameter \code{Asym} is the asymptotic maximum number of species,
  \code{slope} is the maximum slope of increase of richness, and
  \code{xmid} is the  area where half of the maximum richness is
  achieved. 

  In addition to these models, several other models studied by Dengler
  (2009) are available in standard \R self-starting models:
  Michaelis-Menten (\code{\link{SSmicmen}}), Gompertz
  (\code{\link{SSgompertz}}), logistic (\code{\link{SSlogis}}), Weibull
  (\code{\link{SSweibull}}), and some others that may be useful.
    
}
\value{
  Numeric vector of the same length as \code{area}. It is the value of
  the expression of each model. If all arguments are names of objects
  the gradient matrix with respect to these names is attached as an
  attribute named \code{gradient}. 
}
\references{
  Dengler, J. (2009) Which function describes the species-area
  relationship best? A review and empirical evaluation. \emph{Journal of
    Biogeography} 36, 728--744.

  Gitay, H., Roxburgh, S.H. & Wilson, J.B. (1991) Species-area
  relationship in a New Zealand tussock grassland, with implications for
  nature reserve design and for community structure. \emph{Journal of
  Vegetation Science} 2, 113--118.

  Lomolino, M. V. (2000) Ecology's most general, yet protean pattern:
  the species-area relationship. \emph{Journal of Biogeography} 27,
  17--26. 
}
\author{
  Jari Oksanen.
}

\seealso{
  \code{\link{nls}}, \code{\link{fitspecaccum}}. 
}
\examples{
## Get species area data: sipoo.area gives the areas of islands
example(sipoo)
S <- specnumber(sipoo)
plot(S ~ sipoo.area, xlab = "Island Area (ha)", ylab = "Number of Species",
    ylim = c(1, max(S)))
## The Arrhenius model
marr <- nls(S ~ SSarrhenius(sipoo.area, k, z))
marr
## confidence limits from profile likelihood
confint(marr)
## draw a line
xtmp <- seq(min(sipoo.area), max(sipoo.area), len=51)
lines(xtmp, predict(marr, newdata=data.frame(sipoo.area = xtmp)), lwd=2)
## The normal way is to use linear regression on log-log data,
## but this will be different from the previous:
mloglog <- lm(log(S) ~ log(sipoo.area))
mloglog
lines(xtmp, exp(predict(mloglog, newdata=data.frame(sipoo.area=xtmp))),
   lty=2)
## Gleason: log-linear
mgle <- nls(S ~ SSgleason(sipoo.area, k, slope))
lines(xtmp, predict(mgle, newdata=data.frame(sipoo.area=xtmp)),
  lwd=2, col=2)
## Gitay: quadratic of log-linear
mgit <- nls(S ~ SSgitay(sipoo.area, k, slope))
lines(xtmp, predict(mgit, newdata=data.frame(sipoo.area=xtmp)), 
  lwd=2, col = 3)
## Lomolino: using original names of the parameters (Lomolino 2000):
mlom <- nls(S ~ SSlomolino(sipoo.area, Smax, A50, Hill))
mlom
lines(xtmp, predict(mlom, newdata=data.frame(sipoo.area=xtmp)), 
  lwd=2, col = 4)
## One canned model of standard R:
mmic <- nls(S ~ SSmicmen(sipoo.area, slope, Asym))
lines(xtmp, predict(mmic, newdata = data.frame(sipoo.area=xtmp)),
  lwd =2, col = 5)
legend("bottomright", c("Arrhenius", "log-log linear", "Gleason", "Gitay", 
  "Lomolino", "Michaelis-Menten"), col=c(1,1,2,3,4,5), lwd=c(2,1,2,2,2,2), 
   lty=c(1,2,1,1,1,1))
## compare models (AIC)
allmods <- list(Arrhenius = marr, Gleason = mgle, Gitay = mgit, 
   Lomolino = mlom, MicMen= mmic)
sapply(allmods, AIC)
}
\keyword{ models }


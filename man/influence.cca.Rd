\name{influence.cca}
\alias{hatvalues.cca}
\alias{hatvalues.rda}
\alias{sigma.cca}
\alias{rstandard.cca}
\alias{rstudent.cca}
\alias{cooks.distance.cca}
\alias{SSD.cca}
\alias{vcov.cca}
\alias{qr.cca}
\alias{df.residual.cca}

\title{Linear Model Diagnostics for Constrained Ordination}

\description{

  This set of function extracts influence statistics and some other
  linear model statistics directly from a constrained ordination result
  object from \code{\link{cca}}, \code{\link{rda}},
  \code{\link{capscale}} or \code{\link{dbrda}}. The constraints are
  linear model functions and these support functions return identical
  results as the corresponding linear models (\code{\link{lm}}), and you
  can use their documentation. The main functions for normal usage are
  leverage values (\code{\link{hatvalues}}), standardized residuals
  (\code{\link{rstandard}}), studentized or leave-one-out residuals
  (\code{\link{rstudent}}), and Cook's distance
  (\code{\link{cooks.distance}}).  In addition, \code{\link{vcov}}
  returns the variance-covariance matrix of coefficients, and its
  diagonal values the variances of coefficients. Other functions are
  mainly support functions for these, but they can be used directly.

}

\usage{

\method{hatvalues}{cca}(model, ...)
\method{rstandard}{cca}(model, type = c("response", "canoco"), ...)
\method{rstudent}{cca}(model, type = c("response", "canoco"), ...)
\method{cooks.distance}{cca}(model, type = c("response", "canoco"), ...)

\method{sigma}{cca}(object, type = c("response", "canoco"), ...)
\method{vcov}{cca}(object, type = "canoco", ...)
\method{SSD}{cca}(object, type = "canoco", ...)

\method{qr}{cca}(x, ...)
\method{df.residual}{cca}(object, ...)

}

\arguments{

  \item{model, object, x}{A constrained ordination result object.}

  \item{type}{Type of statistics used for extracting raw residuals and
    residual standard deviation (\code{sigma}). Either
    \code{"response"} for species data or difference of WA and LC
    scores for \code{"canoco"}.}

  \item{\dots}{Other arguments to functions (ignored).}

}

\details{

  The \pkg{vegan} algorithm for constrained ordination uses linear model
  (or weighted linear model in \code{\link{cca}}) to find the fitted
  values of dependent community data, and constrained ordination is
  based on this fitted response (Legendre & Legendre 2012). The
  \code{\link{hatvalues}} give the leverage values of these constraints,
  and the leverage is independent on the response data. Other influence
  statistics (\code{\link{rstandard}}, \code{\link{rstudent}},
  \code{\link{cooks.distance}}) are based on leverage, and on the raw
  residuals and residual standard deviation (\code{\link{sigma}}). With
  \code{type = "response"} the raw residuals are given by the
  unconstrained component of the constrained ordination, and influence
  statistics are a matrix with dimensions no. of observations times
  no. of species. For \code{\link{cca}} the statistics are the same as
  obtained from the \code{\link{lm}} model using Chi-square standardized
  species data (see \code{\link{decostand}}) as dependent variable, and
  row sums of community data as weights, and for \code{\link{rda}} the
  \code{\link{lm}} model uses non-modified community data and no
  weights.

  The algorithm in the CANOCO software constraints the results during
  iteration by performing a linear regression of weighted averages (WA)
  scores on constraints and taking the fitted values of this regression
  as linear combination (LC) scores (ter Braak 1984). The WA scores are
  directly found from species scores, but LC scores are linear
  combinations of constraints in the regression. With \code{type =
  "canoco"} the raw residuals are the differences of WA and LC scores,
  and the residual standard deviation (\code{\link{sigma}}) is taken to
  be the axis sum of squared WA scores minus one. These quantities have
  no relationship to residual component of ordination, but they rather
  are methodological artefacts of an algorithm that is not used in
  \pkg{vegan}. The result is a matrix with dimensions no. of
  observations times no. of constrained axes.

  Function \code{\link{vcov}} returns the matrix of variances and
  covariances of regression coefficients. The diagonal values of this
  matrix are the variances, and their square roots give the standard
  errors of regression coefficients. The function is based on
  \code{\link{SSD}} that extracts the sum of squares and crossproducts
  of residuals. The residuals are defined similarly as in influence
  measures and with each \code{type} they have similar properties and
  limitations, and define the dimensions of the result matrix. 

}

\references{

  Legendre, P. and Legendre, L. (2012) \emph{Numerical Ecology}. 3rd
  English ed. Elsevier.
  
  ter Braak, C.J.F. (1984--): CANOCO -- a FORTRAN program for
  \emph{cano}nical \emph{c}ommunity \emph{o}rdination by [partial]
  [detrended] [canonical] correspondence analysis, principal components
  analysis and redundancy analysis. \emph{TNO Inst. of Applied Computer
  Sci., Stat. Dept. Wageningen, The Netherlands}.

  
}

\note{

  Function \code{\link{as.mlm}} casts an ordination object to a multiple
  linear model of class \code{"mlm"} (see \code{\link{lm}}), and similar
  statistics can be derived from that modified object as with this set
  of functions. However, there are some problems in the \R{}
  implementation of the further analysis of multiple linear model
  objects. When the results differ, the current set of functions is more
  probable to be correct. The use of \code{as.mlm} objects should be
  avoided.
  
}

\author{Jari Oksanen}

\seealso{Corresponding \code{\link{lm}} methods and
  \code{\link{as.mlm.cca}}. Function \code{\link{ordiresids}} provides
  lattice graphics for residuals.}

\examples{

data(varespec, varechem)
mod <- cca(varespec ~ Al + P + K, varechem)
## leverage
hatvalues(mod)
plot(hatvalues(mod), type = "h")
## ordination plot with leverages: points with high leverage have
## similar LC and WA scores
plot(mod, type = "n")
ordispider(mod)       # segment from LC to WA scores
points(mod, dis="si", cex=5*hatvalues(mod), pch=21, bg=2) # WA scores
text(mod, dis="bp", col=4)

## deviation and influence
head(rstandard(mod))
head(cooks.distance(mod))

## Influence measures from lm
y <- decostand(varespec, "chi.square") # needed in cca
y1 <- with(y, Cladstel)         # take one species for lm
lmod1 <- lm(y1 ~ Al + P + K, varechem, weights = rowSums(varespec))
## numerically identical within 2e-15
range(cooks.distance(lmod1) - cooks.distance(mod)[, "Cladstel"])

## t-values of regression coefficients based on type = "canoco"
## residuals
coef(mod)
coef(mod)/sqrt(diag(vcov(mod, type = "canoco")))
}

\keyword{ models }
\keyword{ multivariate }
\name{ordisurf}
\alias{ordisurf}
\alias{ordisurf.default}
\alias{ordisurf.formula}
\alias{calibrate.ordisurf}
\alias{plot.ordisurf}

\title{ Fit and Plot Smooth Surfaces of Variables on Ordination. }
\description{
  Function \code{ordisurf} fits a smooth surface for given variable and
  plots the result on ordination diagram.
}
\usage{
\method{ordisurf}{default}(x, y, choices=c(1, 2), knots=10, family="gaussian", col="red",
     thinplate = TRUE, add = FALSE, display = "sites",
     w = weights(x), main, nlevels = 10, levels, labcex = 0.6,
     bubble = FALSE, cex = 1, select = FALSE, method = "GCV.Cp",
     gamma = 1, plot = TRUE, ...)

\method{ordisurf}{formula}(formula, data, ...)

\method{calibrate}{ordisurf}(object, newdata, ...)

\method{plot}{ordisurf}(x, what = c("contour","persp","gam"),
     add = FALSE, bubble = FALSE, col = "red", cex = 1,
     nlevels = 10, levels, labcex = 0.6, \dots)
}

\arguments{
  \item{x}{For \code{ordisurf} an ordination configuration, either a
    matrix or a result known by \code{\link{scores}}. For
    \code{plot.ordisurf} and object of class \code{"ordisurf"} as
    returned by \code{ordisurf}.}
  \item{y}{ Variable to be plotted. }
  \item{choices}{Ordination axes. }
  \item{knots}{Number of initial knots in \code{\link[mgcv]{gam}} (one
    more than degrees of freedom). If \code{knots = 0} or
    \code{knots = 1}  the function will fit a linear trend surface, and
    if \code{knots = 2} the function  will fit a quadratic trend surface
    instead of a smooth surface. } 
  \item{family}{ Error distribution in  \code{\link[mgcv]{gam}}. }
  \item{col}{ Colour of contours. }
  \item{thinplate}{Use thinplate splines in \code{\link[mgcv]{gam}}.}
  \item{add}{Add contours on an existing diagram or draw a new plot. }
  \item{display}{Type of scores known by \code{\link{scores}}: typically
    "sites" for ordinary site scores or "lc" for linear combination scores.}
  \item{w}{Prior weights on the data. Concerns mainly \code{\link{cca}}
    and \code{\link{decorana}} results which have nonconstant weights.}
  \item{main}{The main title for the plot, or as default the name of
    plotted variable in a new plot.}
  \item{nlevels, levels}{Either a vector of \code{levels} for which contours
    are drawn, or suggested number of contours in
    \code{nlevels} if \code{levels} are not supplied.}
  \item{labcex}{Label size in contours.  Setting this zero will suppress
    labels.}
  \item{bubble}{Use \dQuote{bubble plot} for points, or vary the point
    diameter by the value of the plotted variable. If \code{bubble} is
    numeric, its value is used for the maximum symbol size (as in
    \code{cex}), or if \code{bubble = TRUE}, the value of \code{cex} gives
    the maximum. The minimum size will always be \code{cex = 0.4}.  The
    option only has an effect if \code{add = FALSE}.}
  \item{cex}{Character expansion of plotting symbols.}
  \item{select}{Logical; specify \code{\link[mgcv]{gam}} argument
    \code{"select"}. If this is \code{TRUE} then \code{\link[mgcv]{gam}} can
    add an extra  penalty to each term so that it can be penalized to
    zero. This means that the smoothing parameter estimation that is part
    of fitting can completely remove terms from the model. If the
    corresponding smoothing parameter is estimated as zero then the extra
    penalty has no effect.}
  \item{method}{character; the smoothing parameter estimation
    method. Options allowed are: \code{"GCV.Cp"} uses GCV for models with
    unknown scale parameter and Mallows' Cp/UBRE/AIC for models with
    known scale; \code{"GACV.Cp"} as for \code{"GCV.Cp"} but uses GACV
    (Generalised Approximate CV) instead of GCV; \code{"REML"} and
    \code{"ML"} use restricted maximum likelihood or maximum likelihood
    estimation for both known and unknown scale; and \code{"P-REML"} and
    \code{"P-ML"} use REML or ML estimation but use a Pearson estimate
    of the scale.}
  \item{gamma}{Multiplier to inflate model degrees of freedom in GCV or
    UBRE/AIC score by. This effectively places an extra penalty on
    complex models. An oft used value if \code{gamma = 1.4}.}
  \item{plot}{logical; should any plotting be done by
    \code{ordisurf}? Useful if all you want is the fitted response
    surface model.}
  \item{formula, data}{Alternative definition of the fitted model as
    \code{x ~ y}, or left-hand side is the ordination \code{x} and
    right-hand side the single fitted continuous variable
    \code{y}. The variable \code{y} must be in the working environment
    or in the data frame or environment given by \code{data}. All
    other arguments of are passed to the default method.}
  \item{object}{An \code{ordisurf} result object.}
  \item{newdata}{Coordinates in two-dimensional ordination for new
    points.}
  \item{what}{character; what type of plot to produce. \code{"contour"}
    produces a contour plot of the response surface, see
    \code{\link{contour}} for details. \code{"persp"} produces a
    perspective plot of the same, see \code{\link{persp}} for
    details. \code{"gam"} plots the fitted GAM model, an object that
    inherits from class \code{"gam"} returned by \code{ordisurf}, see
    \code{\link[mgcv]{plot.gam}}.} 
  \item{\dots}{Other parameters passed to \code{\link[mgcv]{gam}}, or
    to the graphical functions. See Note below for exceptions.}
}

\details{
  Function \code{ordisurf} fits a smooth surface using thinplate splines
  in \code{\link[mgcv]{gam}}, and uses \code{\link[mgcv]{predict.gam}}
  to find fitted values in a regular grid. The smooth surface can be
  fitted with an extra penalty that allows the entire smoother to be
  penalized back to 0 degrees of freedom, effectively removing the term
  from the model. The addition of this extra penalty is invoked by
  setting argument \code{select} to \code{TRUE}. The function plots the
  fitted contours with convex hull of data points either over an
  existing ordination diagram or draws a new plot. If \code{select ==
    TRUE} and the smooth is effectively penalised out of the model, no
  contours will be plotted.

  \code{\link[mgcv]{gam}} determines the degree of smoothness for the
  fitted response surface during model fitting. Argument \code{method}
  controls how \code{\link[mgcv]{gam}} performs this smoothness
  selection. See \code{\link[mgcv]{gam}} for details of the available
  options. Using \code{"REML"} or \code{"ML"} yields p-values for
  smooths with the best coverage properties if such things matter to
  you.

  The function uses \code{\link{scores}} to extract ordination scores,
  and \code{x} can be any result object known by that function.

  User can supply a vector of prior weights \code{w}. If the ordination
  object has weights, these will be used. In practise this means that
  the row totals are used as weights with
  \code{\link{cca}} or
  \code{\link{decorana}} results. If you do not like this, but want to give
  equal weights to all sites, you should set \code{w = NULL}. The
  behaviour is consistent with \code{\link{envfit}}. For complete
  accordance with constrained \code{\link{cca}}, you should set
  \code{display = "lc"} (and possibly \code{scaling = 2}).

  Function \code{calibrate} returns the fitted values of the response
  variable. The \code{newdata} must be coordinates of points for which
  the fitted values are desired. The function is based on
  \code{\link[mgcv]{predict.gam}} and will pass extra arguments to
  that function.
}

\value{
  Function is usually called for its side effect of drawing the
  contour plot. The function returns the result object of class
  \code{"ordisurf"} that inherits from \code{\link[mgcv]{gam}} used
  internally to fit the surface, but adds an item \code{grid} that
  contains the data for the grid surface. The item \code{grid} has
  elements \code{x} and \code{y} which are vectors of axis coordinates,
  and element \code{z} that is a matrix of fitted values for
  \code{\link{contour}}. The values outside the convex hull of observed
  points are \code{NA} in \code{z}. The \code{\link[mgcv]{gam}}
  component of the result can be used for further analysis like
  predicting new values (see \code{\link[mgcv]{predict.gam}}).
}

\author{ Dave Roberts, Jari Oksanen and Gavin L. Simpson }
\note{   
  The default is to use thinplate splines.  These make sense in
  ordination as they have equal smoothing in all directions and are
  rotation invariant.

  Graphical arguments supplied to \code{plot.ordisurf} are passed on to
  the underlying plotting functions, \code{contour}, \code{persp}, and
  \code{\link[mgcv]{plot.gam}}. The exception to this is that arguments
  \code{col} and \code{cex} can not currently be passed to
  \code{\link[mgcv]{plot.gam}} because of a bug in the way that function
  evaluates arguments when arranging the plot.

  A work-around is to call \code{\link[mgcv]{plot.gam}} directly on the
  result of a call to \code{ordisurf}. See the Examples for an
  illustration of this.  
}

\seealso{ For basic routines \code{\link[mgcv]{gam}},
  and \code{\link{scores}}. Function 
  \code{\link{envfit}} provides a more traditional and compact
  alternative. } 

\examples{
data(varespec)
data(varechem)
library(MASS)
vare.dist <- vegdist(varespec)
vare.mds <- isoMDS(vare.dist)
with(varechem, ordisurf(vare.mds, Baresoil, bubble = 5))

## as above but with extra penalties on smooth terms:
with(varechem, ordisurf(vare.mds, Baresoil, bubble = 5, col = "blue",
                        add = TRUE, select = TRUE))

## Cover of Cladina arbuscula
fit <- with(varespec, ordisurf(vare.mds, Cla.arb, family=quasipoisson)) 
## Get fitted values
calibrate(fit)

## Plot method
plot(fit, what = "contour")

## Plotting the "gam" object
plot(fit, what = "gam") ## 'col' and 'cex' not passed on
## or via plot.gam directly
plot.gam(fit, cex = 2, pch = 1, col = "blue")
## 'col' effects all objects drawn...
}
\keyword{ multivariate }
\keyword{ aplot }

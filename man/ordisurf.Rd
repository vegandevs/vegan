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
\method{ordisurf}{default}(x, y, choices = c(1, 2), knots = 10,
         family = "gaussian", col = "red", isotropic = TRUE,
         thinplate = TRUE, bs = "tp", fx = FALSE, add = FALSE,
         display = "sites", w = weights(x, display), main, nlevels = 10,
         levels, npoints = 31, labcex = 0.6, bubble = FALSE,
         cex = 1, select = TRUE, method = "REML", gamma = 1,
         plot = TRUE, lwd.cl = par("lwd"), ...)

\method{ordisurf}{formula}(formula, data, ...)

\method{calibrate}{ordisurf}(object, newdata, ...)

\method{plot}{ordisurf}(x, what = c("contour","persp","gam"),
     add = FALSE, bubble = FALSE, col = "red", cex = 1,
     nlevels = 10, levels, labcex = 0.6, lwd.cl = par("lwd"), \dots)
}

\arguments{
  \item{x}{For \code{ordisurf} an ordination configuration, either a
    matrix or a result known by \code{\link{scores}}. For
    \code{plot.ordisurf} an object of class \code{"ordisurf"} as
    returned by \code{ordisurf}.}
  \item{y}{Variable to be plotted / modelled as a function of the
    ordination scores.}
  \item{choices}{Ordination axes. }
  \item{knots}{Number of initial knots in \code{\link[mgcv]{gam}} (one
    more than degrees of freedom). If \code{knots = 0} or
    \code{knots = 1}  the function will fit a linear trend surface, and
    if \code{knots = 2} the function  will fit a quadratic trend surface
    instead of a smooth surface. A vector of length 2 is allowed when
    \code{isotropic = FALSE}, with the first and second elements of
    \code{knots} referring to the first and second of ordination
    dimensions (as indicated by \code{choices}) respectively.}
  \item{family}{Error distribution in \code{\link[mgcv]{gam}}.}
  \item{col}{ Colour of contours. }
  \item{isotropic, thinplate}{Fit an isotropic smooth surface (i.e. same
    smoothness in both ordination dimensions) via
    \code{\link[mgcv]{gam}}. Use of \code{thinplate} is deprecated and
    will be removed in a future version of the package.}
  \item{bs}{a two letter character string indicating the smoothing basis
    to use. (e.g. \code{"tp"} for thin plate regression spline,
    \code{"cr"} for cubic regression spline). One of \code{c("tp", "ts",
      "cr", "cs", "ds", "ps", "ad")}. See
    \code{\link[mgcv]{smooth.terms}} for an over view of what these
    refer to. The default is to use thin plate splines: \code{bs = "tp"}.}
  \item{fx}{indicates whether the smoothers are fixed degree of freedom
    regression splines (\code{fx = FALSE}) or penalised regression
    splines (\code{fx = TRUE}). Can be a vector of length 2 for
    anisotropic surfaces (\code{isotropic = FALSE}). It doesn't make
    sense to use \code{fx = TRUE} \strong{and} \code{select = TRUE} and
    it is an \strong{error} to do so. A warning is issued if you specify
    \code{fx = TRUE} and forget to use \code{select = FALSE} though
    fitting continues using \code{select = FALSE}.}
  \item{add}{Add contours to an existing diagram or draw a new plot?}
  \item{display}{Type of scores known by \code{\link{scores}}: typically
    "sites" for ordinary site scores or "lc" for linear combination scores.}
  \item{w}{Prior weights on the data. Concerns mainly \code{\link{cca}}
    and \code{\link{decorana}} results which have nonconstant weights.}
  \item{main}{The main title for the plot, or as default the name of
    plotted variable in a new plot.}
  \item{nlevels, levels}{Either a vector of \code{levels} for which contours
    are drawn, or suggested number of contours in \code{nlevels} if
    \code{levels} are not supplied.}
  \item{npoints}{numeric; the number of locations at which to evaluate
    the fitted surface. This represents the number of locations in each
    dimension.}
  \item{labcex}{Label size in contours.  Setting this zero will suppress
    labels.}
  \item{bubble}{Use a \dQuote{bubble plot} for points, or vary the point
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
    complex models. An oft-used value is \code{gamma = 1.4}.}
  \item{plot}{logical; should any plotting be done by
    \code{ordisurf}? Useful if all you want is the fitted response
    surface model.}
  \item{lwd.cl}{numeric; the \code{lwd} (line width) parameter to use
    when drawing the contour lines.}
  \item{formula, data}{Alternative definition of the fitted model as
    \code{x ~ y}, where left-hand side is the ordination \code{x} and
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
  \item{\dots}{Other parameters passed to \code{\link{scores}}, or
    to the graphical functions. See Note below for exceptions.}
}

\details{

  Function \code{ordisurf} fits a smooth surface using penalised
  splines (Wood 2003) in \code{\link[mgcv]{gam}}, and uses
  \code{\link[mgcv]{predict.gam}} to find fitted values in a regular
  grid. The smooth surface can be fitted with an extra penalty that
  allows the entire smoother to be penalized back to 0 degrees of
  freedom, effectively removing the term from the model (see Marra &
  Wood, 2011). The addition of this extra penalty is invoked by
  setting argument \code{select} to \code{TRUE}. An alternative is to
  use a spline basis that includes shrinkage (\code{bs = "ts"} or
  \code{bs = "cs"}).

  \code{ordisurf()} exposes a large number of options from
  \code{\link[mgcv]{gam}} for specifying the basis functions used for
  the surface. If you stray from the defaults, do read the
  \strong{Notes} section below and relevant documentation in
  \code{\link[mgcv]{s}} and \code{\link[mgcv]{smooth.terms}}.

  The function plots the fitted contours with convex hull of data points
  either over an existing ordination diagram or draws a new plot. If
  \code{select = TRUE} and the smooth is effectively penalised out of
  the model, no contours will be plotted.

  \code{\link[mgcv]{gam}} determines the degree of smoothness for the
  fitted response surface during model fitting, unless \code{fx =
  TRUE}. Argument \code{method} controls how \code{\link[mgcv]{gam}}
  performs this smoothness selection. See \code{\link[mgcv]{gam}} for
  details of the available options. Using \code{"REML"} or \code{"ML"}
  yields p-values for smooths with the best coverage properties if such
  things matter to you.

  The function uses \code{\link{scores}} to extract ordination scores,
  and \code{x} can be any result object known by that function.

  The user can supply a vector of prior weights \code{w}. If the
  ordination object has weights, these will be used. In practise this
  means that the row totals are used as weights with \code{\link{cca}}
  or \code{\link{decorana}} results. If you do not like this, but want
  to give equal weights to all sites, you should set \code{w =
  NULL}. The behaviour is consistent with \code{\link{envfit}}. For
  complete accordance with constrained \code{\link{cca}}, you should set
  \code{display = "lc"}.

  Function \code{calibrate} returns the fitted values of the response
  variable. The \code{newdata} must be coordinates of points for which
  the fitted values are desired. The function is based on
  \code{\link[mgcv]{predict.gam}} and will pass extra arguments to
  that function.
}

\value{
  \code{ordisurf} is usually called for its side effect of drawing the
  contour plot. The function returns a result object of class
  \code{"ordisurf"} that inherits from \code{\link[mgcv]{gam}} used
  internally to fit the surface, but adds an item \code{grid} that
  contains the data for the grid surface. The item \code{grid} has
  elements \code{x} and \code{y} which are vectors of axis coordinates,
  and element \code{z} that is a matrix of fitted values for
  \code{\link{contour}}. The values outside the convex hull of observed
  points are indicated as \code{NA} in \code{z}. The
  \code{\link[mgcv]{gam}} component of the result can be used for
  further analysis like predicting new values (see
  \code{\link[mgcv]{predict.gam}}).
}

\author{ Dave Roberts, Jari Oksanen and Gavin L. Simpson }
\note{
  The default is to use an isotropic smoother via
  \code{\link[mgcv]{s}} employing thin plate regression splines
  (\code{bs = "tp"}). These make sense in ordination as they have
  equal smoothing in all directions and are rotation invariant. However,
  if different degrees of smoothness along dimensions are required, an
  anisotropic smooth surface may be more applicable. This can be
  achieved through the use of \code{isotropic = FALSE}, wherein the
  surface is fitted via a tensor product smoother via
  \code{\link[mgcv]{te}} (unless \code{bs = "ad"}, in which case
  separate splines for each dimension are fitted using
  \code{\link[mgcv]{s}}).

  Cubic regression splines and P splines can \strong{only} be used with
  \code{isotropic = FALSE}.

  Adaptive smooths (\code{bs = "ad"}), especially in two dimensions,
  require a large number of observations; without many hundreds of
  observations, the default complexities for the smoother will exceed
  the number of observations and fitting will fail.

  To get the old behaviour of \code{ordisurf} use \code{select = FALSE},
  \code{method = "GCV.Cp"}, \code{fx = FALSE}, and \code{bs = "tp"}. The
  latter two options are the current defaults.

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

\section{Warning}{
  The fitted GAM is a regression model and has the usual assumptions of
  such models. Of particular note is the assumption of independence of
  residuals. If the observations are not independent (e.g. they are
  repeat measures on a set of objects, or from an experimental design,
  \emph{inter alia}) do not trust the \emph{p}-values from the GAM
  output.

  If you need further control (i.e. to add additional fixed effects to
  the model, or use more complex smoothers), extract the ordination
  scores using the \code{scores} function and then generate your own
  \code{\link[mgcv]{gam}} call.
}

\references{

  Marra, G.P & Wood, S.N. (2011) Practical variable selection for
  generalized additive models. \emph{Comput. Stat. Data Analysis} 55,
  2372--2387.

  Wood, S.N. (2003) Thin plate regression splines.
  \emph{J. R. Statist. Soc. B} 65, 95--114.

}

\seealso{ For basic routines \code{\link[mgcv]{gam}},
  and \code{\link{scores}}. Function
  \code{\link{envfit}} provides a more traditional and compact
  alternative. }

\examples{
data(varespec)
data(varechem)
vare.dist <- vegdist(varespec)
vare.mds <- monoMDS(vare.dist)
ordisurf(vare.mds ~ Baresoil, varechem, bubble = 5)

## as above but without the extra penalties on smooth terms,
## and using GCV smoothness selection (old behaviour of `ordisurf()`):
ordisurf(vare.mds ~ Baresoil, varechem, col = "blue", add = TRUE,
                        select = FALSE, method = "GCV.Cp")

## Cover of Cladina arbuscula
fit <- ordisurf(vare.mds ~ Cladarbu, varespec, family=quasipoisson)
## Get fitted values
calibrate(fit)

## Variable selection via additional shrinkage penalties
## This allows non-significant smooths to be selected out
## of the model not just to a linear surface. There are 2
## options available:
##  - option 1: `select = TRUE` --- the *default*
ordisurf(vare.mds ~ Baresoil, varechem, method = "REML", select = TRUE)
##  - option 2: use a basis with shrinkage
ordisurf(vare.mds ~ Baresoil, varechem, method = "REML", bs = "ts")
## or bs = "cs" with `isotropic = FALSE`

## Plot method
plot(fit, what = "contour")

## Plotting the "gam" object
plot(fit, what = "gam") ## 'col' and 'cex' not passed on
## or via plot.gam directly
library(mgcv)
plot.gam(fit, cex = 2, pch = 1, col = "blue")
## 'col' effects all objects drawn...

### controlling the basis functions used
## Use Duchon splines
ordisurf(vare.mds ~ Baresoil, varechem, bs = "ds")

## A fixed degrees of freedom smooth, must use 'select = FALSE'
ordisurf(vare.mds ~ Baresoil, varechem, knots = 4,
                        fx = TRUE, select = FALSE)

## An anisotropic smoother with cubic regression spline bases
ordisurf(vare.mds ~ Baresoil, varechem, isotropic = FALSE,
                        bs = "cr", knots = 4)

## An anisotropic smoother with cubic regression spline with
## shrinkage bases & different degrees of freedom in each dimension
ordisurf(vare.mds ~ Baresoil, varechem, isotropic = FALSE,
                        bs = "cs", knots = c(3,4), fx = TRUE,
                        select = FALSE)
}
\keyword{ multivariate }
\keyword{ aplot }

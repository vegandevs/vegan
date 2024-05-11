\name{decorana}
\alias{decorana}
\alias{summary.decorana}
\alias{print.summary.decorana}
\alias{plot.decorana}
\alias{downweight}
\alias{scores.decorana}
\alias{points.decorana}
\alias{text.decorana}

\title{Detrended Correspondence Analysis and Basic Reciprocal Averaging }
\description{
  Performs detrended correspondence analysis and basic reciprocal
  averaging or orthogonal correspondence analysis.
}
\usage{
decorana(veg, iweigh=0, iresc=4, ira=0, mk=26, short=0,
         before=NULL, after=NULL)

\method{plot}{decorana}(x, choices=c(1,2), origin=TRUE,
     display=c("both","sites","species","none"),
     cex = 0.8, cols = c(1,2), type, xlim, ylim, ...)

\method{text}{decorana}(x, display = c("sites", "species"), labels,
     choices = 1:2, origin = TRUE, select,  ...)

\method{points}{decorana}(x, display = c("sites", "species"),
       choices=1:2, origin = TRUE, select, ...)

## deprecated: use 'scores' for scores, 'weights' for weights
\method{summary}{decorana}(object, digits=3, origin=TRUE,
        display=c("both", "species","sites","none"), ...)

\method{print}{summary.decorana}(x, head = NA, tail = head, ...)

downweight(veg, fraction = 5)

\method{scores}{decorana}(x, display="sites", choices=1:4,
       origin=TRUE, tidy=FALSE, ...)
}

\arguments{
  \item{veg}{Community data, a matrix-like object. }
  \item{iweigh}{Downweighting of rare species (0: no). }
  \item{iresc}{Number of rescaling cycles (0: no rescaling). }
  \item{ira}{Type of analysis (0: detrended, 1: basic reciprocal averaging). }
  \item{mk}{Number of segments in rescaling. }
  \item{short}{Shortest gradient to be rescaled. }
  \item{before}{Hill's piecewise transformation: values before transformation. }
  \item{after}{Hill's piecewise transformation: values after
  transformation -- these must correspond to values in \code{before}.}
  \item{x, object}{A \code{decorana} result object.}
  \item{choices}{Axes shown.}
  \item{origin}{Use true origin even in detrended correspondence analysis.}
  \item{display}{Display only sites, only species, both or neither.}
  \item{cex}{Plot character size.}
  \item{cols}{Colours used for sites and species.}
  \item{type}{Type of plots, partial match to \code{"text"},
    \code{"points"} or \code{"none"}.}
  \item{labels}{Optional text to be used instead of row names.}
  \item{select}{Items to be displayed.  This can either be a logical
    vector which is \code{TRUE} for displayed items or a vector of indices
    of displayed items.}
  \item{xlim, ylim}{the x and y limits (min,max) of the plot.}
  \item{digits}{Number of digits in summary output.}
  \item{head, tail}{Number of rows printed from the head and tail of
    species and site scores. Default \code{NA} prints all.}
  \item{fraction}{Abundance fraction where downweighting begins.}
  \item{tidy}{Return scores that are compatible with \CRANpkg{ggplot2}:
    all scores are in a single \code{data.frame}, score type is
    identified by factor variable \code{score} (\code{"sites"},
    \code{"species"}), the names by variable \code{label}. These scores
    are incompatible with conventional \code{plot} functions, but they can
    be used in \pkg{ggplot2}.}
  \item{...}{Other arguments for \code{plot} function.}
  }

\details{
  In late 1970s, correspondence analysis became the method of choice for
  ordination in vegetation science, since it seemed better able to cope 
  with non-linear species responses than principal components
  analysis. However, even correspondence analysis can produce an arc-shaped
  configuration of a single gradient. Mark Hill developed detrended
  correspondence analysis to correct two assumed \sQuote{faults} in 
  correspondence analysis: curvature of straight gradients and packing
  of sites at the ends of the gradient.  

  The curvature is removed by replacing the orthogonalization of axes
  with detrending.  In orthogonalization successive axes are made
  non-correlated, but detrending should remove all systematic dependence
  between axes.  Detrending is performed using a smoothing window on
  \code{mk} segments.  The packing of sites at the ends of the gradient
  is undone by rescaling the axes after extraction.  After rescaling,
  the axis is supposed to be scaled by \sQuote{SD} units, so that the
  average width of Gaussian species responses is supposed to be one over
  whole axis. Other innovations were the piecewise linear transformation
  of species abundances and downweighting of rare species which were
  regarded to have an unduly high influence on ordination axes.

  It seems that detrending actually works by twisting the ordination
  space, so that the results look non-curved in two-dimensional
  projections (\sQuote{lolly paper effect}).  As a result, the points
  usually have an easily recognized triangular or diamond shaped
  pattern, obviously an artefact of detrending.  Rescaling works
  differently than commonly presented, too. \code{decorana} does not
  use, or even evaluate, the widths of species responses.  Instead, it
  tries to equalize the weighted standard deviation of species scores on
  axis segments (parameter \code{mk} has no effect, since
  \code{decorana} finds the segments internally). Function
  \code{\link{tolerance}} returns this internal criterion and can be
  used to assess the success of rescaling.

  The \code{plot} method plots species and site scores.  Classical
  \code{decorana} scaled the axes so that smallest site score was 0 (and
  smallest species score was negative), but \code{summary}, \code{plot}
  and \code{scores} use the true origin, unless \code{origin = FALSE}.

  In addition to proper eigenvalues, the function reports
  \sQuote{decorana values} in detrended analysis. These \sQuote{decorana
  values} are the values that the legacy code of \code{decorana} returns
  as eigenvalues. They are estimated during iteration, and describe the
  joint effects of axes and detrending. The \sQuote{decorana values} are
  estimated before rescaling and do not show its effect on
  eigenvalues. The proper eigenvalues are estimated after extraction of
  the axes and they are the ratio of weighted sum of squares of site and
  species scores even in detrended and rescaled solutions. These
  eigenvalues are estimated for each axis separately, but they are not
  additive, because higher \code{decorana} axes can show effects already
  explained by prior axes. \sQuote{Additive eigenvalues} are cleansed
  from the effects of prior axes, and they can be assumed to add up to
  total inertia (scaled Chi-square). For proportions and cumulative
  proportions explained you can use \code{\link{eigenvals.decorana}}.

}
\value{
  \code{decorana} returns an object of class \code{"decorana"}, which
  has \code{print}, \code{summary}, \code{scores}, \code{plot},
  \code{points} and \code{text} methods, and support functions
  \code{\link{eigenvals}}, \code{\link{bstick}}, \code{screeplot},
  \code{predict} and \code{\link{tolerance}}. \code{downweight} is an
  independent function that can also be used with other methods than
  \code{decorana}.
}

\references{
  Hill, M.O. and Gauch, H.G. (1980). Detrended correspondence analysis:
  an improved ordination technique. \emph{Vegetatio} \strong{42},
  47--58.

  Oksanen, J. and Minchin, P.R. (1997). Instability of ordination
  results under changes in input data order: explanations and
  remedies. \emph{Journal of Vegetation Science} \strong{8}, 447--454.
}
\author{Mark O. Hill wrote the original Fortran code, the \R port was by
  Jari Oksanen. }
\note{
  \code{decorana} uses the central numerical engine of the
  original Fortran code (which is in the public domain), or about 1/3 of
  the original program.  I have tried to implement the original
  behaviour, although a great part of preparatory steps were written in
  \R language, and may differ somewhat from the original code. However,
  well-known bugs are corrected and strict criteria used (Oksanen &
  Minchin 1997). 

  Please note that there really is no need for piecewise transformation
  or even downweighting within \code{decorana}, since there are more
  powerful and extensive alternatives in \R, but these options are
  included for compliance with the original software.  If a different
  fraction of abundance is needed in downweighting, function
  \code{downweight} must be applied before \code{decorana}.  Function
  \code{downweight} indeed can be applied prior to correspondence
  analysis, and so it can be used together with \code{\link{cca}}, too.

  Github package \pkg{natto} has an \R{} implementation of
  \code{decorana} which allows easier inspection of the
  algorithm and also easier development of the function.
}

\seealso{

  For unconstrained ordination, non-metric multidimensional scaling in
  \code{\link{monoMDS}} may be more robust (see also
  \code{\link{metaMDS}}).  Constrained (or \sQuote{canonical})
  correspondence analysis can be made with \code{\link{cca}}.
  Orthogonal correspondence analysis can be made with \code{decorana} or
  \code{\link{cca}}, but the scaling of results vary (and the one in
  \code{decorana} corresponds to \code{scaling = "sites"} and \code{hill
  = TRUE} in \code{\link{cca}}.).  See \code{\link{predict.decorana}}
  for adding new points to an ordination.

}

\examples{
data(varespec)
vare.dca <- decorana(varespec)
vare.dca
plot(vare.dca)

### the detrending rationale:
gaussresp <- function(x,u) exp(-(x-u)^2/2)
x <- seq(0,6,length=15) ## The gradient
u <- seq(-2,8,len=23)   ## The optima
pack <- outer(x,u,gaussresp)
matplot(x, pack, type="l", main="Species packing")
opar <- par(mfrow=c(2,2))
plot(scores(prcomp(pack)), asp=1, type="b", main="PCA")
plot(scores(decorana(pack, ira=1)), asp=1, type="b", main="CA")
plot(scores(decorana(pack)), asp=1, type="b", main="DCA")
plot(scores(cca(pack ~ x), dis="sites"), asp=1, type="b", main="CCA")

### Let's add some noise:
noisy <- (0.5 + runif(length(pack)))*pack
par(mfrow=c(2,1))
matplot(x, pack, type="l", main="Ideal model")
matplot(x, noisy, type="l", main="Noisy model")
par(mfrow=c(2,2))
plot(scores(prcomp(noisy)), type="b", main="PCA", asp=1)
plot(scores(decorana(noisy, ira=1)), type="b", main="CA", asp=1)
plot(scores(decorana(noisy)), type="b", main="DCA", asp=1)
plot(scores(cca(noisy ~ x), dis="sites"), asp=1, type="b", main="CCA")
par(opar)
}
\keyword{ multivariate }


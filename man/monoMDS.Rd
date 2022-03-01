\name{monoMDS}
\alias{monoMDS}
\alias{scores.monoMDS}
\alias{plot.monoMDS}
\alias{points.monoMDS}
\alias{text.monoMDS}

\title{ Global and Local Non-metric Multidimensional Scaling and
  Linear and Hybrid Scaling }

\description{ Function implements Kruskal's (1964a,b) non-metric
  multidimensional scaling (NMDS) using monotone regression and
  primary (\dQuote{weak}) treatment of ties. In addition to
  traditional global NMDS, the function implements local NMDS, linear
  and hybrid multidimensional scaling.
}

\usage{
monoMDS(dist, y, k = 2, model = c("global", "local", "linear", "hybrid"),
    threshold = 0.8, maxit = 200, weakties = TRUE, stress = 1,
    scaling = TRUE, pc = TRUE, smin = 1e-4, sfgrmin = 1e-7,
    sratmax=0.999999, ...)
\method{scores}{monoMDS}(x, choices = NA, ...)
\method{plot}{monoMDS}(x, choices = c(1,2), type = "t", ...)
\method{points}{monoMDS}(x, choices = c(1,2), select, ...)
\method{text}{monoMDS}(x, labels, choices = c(1,2), select, ...)
}

\arguments{
  \item{dist}{Input dissimilarities.}
  \item{y}{Starting configuration. A random configuration will be
    generated if this is missing.}

  \item{k}{Number of dimensions. NB., the number of points \eqn{n}
    should be \eqn{n > 2k + 1}{n > 2*k + 1}, and preferably higher in
    non-metric MDS.}

  \item{model}{MDS model: \code{"global"} is normal non-metric MDS
    with a monotone regression, \code{"local"} is non-metric MDS with
    separate regressions for each point, \code{"linear"} uses linear
    regression, and \code{"hybrid"} uses linear regression for
    dissimilarities below a threshold in addition to monotone
    regression. See Details.}

  \item{threshold}{Dissimilarity below which linear regression is
    used alternately with monotone regression. }

  \item{maxit}{Maximum number of iterations.}

  \item{weakties}{Use primary or weak tie treatment, where equal
    observed dissimilarities are allowed to have different fitted
    values. if \code{FALSE}, then secondary (strong) tie treatment is
    used, and tied values are not broken.}

  \item{stress}{Use stress type 1 or 2 (see Details).}

  \item{scaling}{Scale final scores to unit root mean squares.}

  \item{pc}{Rotate final scores to principal components.}

  \item{smin, sfgrmin, sratmax}{Convergence criteria: iterations stop
    when stress drops below \code{smin}, scale factor of the gradient
    drops below \code{sfgrmin}, or stress ratio between two iterations
    goes over \code{sratmax} (but is still \eqn{< 1}).}

  \item{x}{A \code{monoMDS} result.}

  \item{choices}{Dimensions returned or plotted. The default \code{NA}
    returns all dimensions. }

  \item{type}{The type of the plot: \code{"t"} for text, \code{"p"}
    for points, and \code{"n"} for none.}

  \item{select}{Items to be displayed.  This can either be a logical
          vector which is \code{TRUE} for displayed items or a vector of
          indices of displayed items.}

  \item{labels}{Labels to be use used instead of row names.}

  \item{\dots}{Other parameters to the functions (ignored in
    \code{monoMDS}, passed to graphical functions in \code{plot}.).}
}

\details{ There are several versions of non-metric multidimensional
  scaling in \R, but \code{monoMDS} offers the following unique
  combination of features:
  \itemize{

    \item \dQuote{Weak} treatment of ties (Kruskal 1964a,b), where tied
    dissimilarities can be broken in monotone regression. This is
    especially important for cases where compared sites share no species
    and dissimilarities are tied to their maximum value of one. Breaking
    ties allows these points to be at different distances and can help
    in recovering very long coenoclines (gradients).  Functions in the
    \CRANpkg{smacof} package also hav adequate tie treatment.

    \item Handles missing values in a meaningful way.

    \item Offers \dQuote{local} and \dQuote{hybrid} scaling in
    addition to usual \dQuote{global} NMDS (see below).

    \item Uses fast compiled code (\code{\link[MASS]{isoMDS}} of the
    \pkg{MASS} package also uses compiled code).

 }

  Function \code{monoMDS} uses Kruskal's (1964b) original monotone
  regression to minimize the stress. There are two alternatives of
  stress: Kruskal's (1964a,b) original or \dQuote{stress 1} and an
  alternative version or \dQuote{stress 2} (Sibson 1972). Both of
  these stresses can be expressed with a general formula

\deqn{s^2 = \frac{\sum (d - \hat d)^2}{\sum(d - d_0)^2}}{stress^2 = sum (d-dhat)^2/ sum (d-dnull)^2}

  where \eqn{d} are distances among points in ordination configuration,
  \eqn{\hat d}{dhat} are the fitted ordination distances, and
  \eqn{d_0}{dnull} are the ordination distances under null model.  For
  \dQuote{stress 1} \eqn{d_0 = 0}{dnull = 0}, and for \dQuote{stress 2}
  \eqn{d_0 = \bar{d}}{dnull = dbar} or mean distances. \dQuote{Stress 2}
  can be expressed as \eqn{s^2 = 1 - R^2}{stress^2 = 1 - R2},
  where\eqn{R^2}{R2} is squared correlation between fitted values and
  ordination distances, and so related to the \dQuote{linear fit} of
  \code{\link{stressplot}}.

  Function \code{monoMDS} can fit several alternative NMDS variants that
  can be selected with argument \code{model}.  The default \code{model =
  "global"} fits global NMDS, or Kruskal's (1964a,b) original NMDS
  similar to \code{\link[MASS]{isoMDS}} (\pkg{MASS}).  Alternative
  \code{model = "local"} fits local NMDS where independent monotone
  regression is used for each point (Sibson 1972).  Alternative
  \code{model = "linear"} fits a linear MDS. This fits a linear
  regression instead of monotone, and is not identical to metric scaling
  or principal coordinates analysis (\code{\link{cmdscale}}) that
  performs an eigenvector decomposition of dissimilarities (Gower
  1966). Alternative \code{model = "hybrid"} implements hybrid MDS that
  uses monotone regression for all points and linear regression for
  dissimilarities below or at a \code{threshold} dissimilarity in
  alternating steps (Faith et al. 1987). Function
  \code{\link{stressplot}} can be used to display the kind of regression
  in each \code{model}.

  Scaling, orientation and direction of the axes is arbitrary.
  However, the function always centres the axes, and the default
  \code{scaling} is to scale the configuration of unit root mean
  square and to rotate the axes (argument \code{pc}) to principal
  components so that the first dimension shows the major variation.
  It is possible to rotate the solution so that the first axis is
  parallel to a given environmental variable using function
  \code{\link{MDSrotate}}.

}

\section{Convergence Criteria}{

  NMDS is iterative, and the function stops when any of its
  convergence criteria is met. There is actually no criterion of
  assured convergence, and any solution can be a local optimum. You
  should compare several random starts (or use \code{monoMDS} via
  \code{\link{metaMDS}}) to assess if the solutions is likely a global
  optimum.

  The stopping criteria are:
  \describe{

     \item{\code{maxit}:}{ Maximum number of iterations. Reaching this
     criterion means that solutions was almost certainly not found,
     and \code{maxit} should be increased.}

     \item{\code{smin}:}{ Minimum stress. If stress is nearly zero,
     the fit is almost perfect. Usually this means that data set is
     too small for the requested analysis, and there may be several
     different solutions that are almost as perfect. You should reduce
     the number of dimensions (\code{k}), get more data (more
     observations) or use some other method, such as metric scaling
     (\code{\link{cmdscale}}, \code{\link{wcmdscale}}).}

     \item{\code{sratmax}:}{ Change in stress. Values close to one
     mean almost unchanged stress. This may mean a solution, but it
     can also signal stranding on suboptimal solution with flat stress
     surface.}

     \item{\code{sfgrmin}:}{ Minimum scale factor. Values close to
     zero mean almost unchanged configuration. This may mean a
     solution, but will also happen in local optima.}

  }

}

\value{ Returns an object of class \code{"monoMDS"}. The final scores
  are returned in item \code{points} (function \code{scores} extracts
  these results), and the stress in item \code{stress}. In addition,
  there is a large number of other items (but these may change without
  notice in the future releases). }

\references{

  Faith, D.P., Minchin, P.R and Belbin, L. 1987. Compositional
      dissimilarity as a robust measure of ecological
      distance. \emph{Vegetatio} 69, 57--68.

  Gower, J.C. (1966). Some distance properties of latent root and
      vector methods used in multivariate analysis. \emph{Biometrika}
      53, 325--328.

  Kruskal, J.B. 1964a. Multidimensional scaling by optimizing
      goodness-of-fit to a nonmetric hypothesis. \emph{Psychometrika}
      29, 1--28.

  Kruskal, J.B. 1964b. Nonmetric multidimensional scaling: a numerical
      method. \emph{Psychometrika} 29, 115--129.

  Minchin, P.R. 1987. An evaluation of relative robustness of
      techniques for ecological ordinations. \emph{Vegetatio} 69,
      89--107.

  Sibson, R. 1972. Order invariant methods for data
      analysis. \emph{Journal of the Royal Statistical Society B} 34,
      311--349.

}

\author{
Peter R. Michin (Fortran core) and Jari Oksanen (R interface).
}

\note{ This is the default NMDS function used in
   \code{\link{metaMDS}}. Function \code{\link{metaMDS}} adds support
   functions so that NMDS can be run like recommended by Minchin
   (1987).
}

\seealso{ \code{\link{metaMDS}} for the \pkg{vegan} way of running NMDS,
  and \code{\link[MASS]{isoMDS}} and \CRANpkg{smacof} for some
  alternative implementations of NMDS. }

\examples{
data(dune)
dis <- vegdist(dune)
m <- monoMDS(dis, model = "loc")
m
plot(m)
}
\keyword{ multivariate }


\name{radfit}
\alias{radfit}
\alias{radfit.default}
\alias{radfit.data.frame}
\alias{AIC.radfit}
\alias{AIC.radfit.frame}
\alias{as.rad}
\alias{coef.radfit}
\alias{coef.radfit.frame}
\alias{deviance.radfit}
\alias{deviance.radfit.frame}
\alias{logLik, radfit}
\alias{logLik, radfit.frame}
\alias{fitted.radfit}
\alias{fitted.radfit.frame}
\alias{lines.radline}
\alias{lines.radfit}
\alias{plot.radfit.frame}
\alias{plot.radfit}
\alias{plot.radline}
\alias{plot.rad}
\alias{radlattice}
\alias{points.radline}
\alias{points.radfit}
\alias{summary.radfit.frame}
\alias{rad.preempt}
\alias{rad.lognormal}
\alias{rad.zipf}
\alias{rad.zipfbrot}
\alias{rad.null}
\alias{predict.radline}
\alias{predict.radfit}
\alias{predict.radfit.frame}

\title{ Rank -- Abundance or Dominance / Diversity Models}
\description{
  Functions construct rank -- abundance or dominance / diversity or
  Whittaker plots and fit brokenstick, pre-emption, log-Normal,
  Zipf and Zipf-Mandelbrot models of species abundance.
}
\usage{
\method{radfit}{default}(x, ...)
rad.null(x, family=poisson, ...)
rad.preempt(x, family = poisson, ...)
rad.lognormal(x, family = poisson, ...)
rad.zipf(x, family = poisson, ...)
rad.zipfbrot(x, family = poisson, ...)
\method{predict}{radline}(object, newdata, total, ...)
\method{plot}{radfit}(x, BIC = FALSE, legend = TRUE, ...)  
\method{plot}{radfit.frame}(x, order.by, BIC = FALSE, model, legend = TRUE,
     as.table = TRUE, ...)
\method{plot}{radline}(x, xlab = "Rank", ylab = "Abundance", type = "b", ...)
radlattice(x, BIC = FALSE, ...)
\method{lines}{radfit}(x, ...)
\method{points}{radfit}(x, ...)
as.rad(x)
\method{plot}{rad}(x, xlab = "Rank", ylab = "Abundance", log = "y", ...)
}

\arguments{
  \item{x}{Data frame, matrix or a vector giving species abundances, or an object to
    be plotted.}

  \item{family}{Error distribution (passed to \code{\link{glm}}). All
    alternatives accepting \code{link = "log"} in \code{\link{family}}
    can be used, although not all make sense.}

  \item{object}{A fitted result object.}

  \item{newdata}{Ranks used for ordinations. All models can
    interpolate to non-integer \dQuote{ranks} (although this may be
    approximate), but extrapolation may fail}

  \item{total}{The new total used for predicting abundance. Observed
    total count is used if this is omitted.}

  \item{order.by}{A vector used for ordering sites in plots.}
  \item{BIC}{Use Bayesian Information Criterion, BIC, instead of
    Akaike's AIC. The penalty in BIC is \eqn{k = \log(S)}{k =
      log(S)}  where \eqn{S} is the number of species, whereas AIC uses
    \eqn{k = 2}.} 
  \item{model}{Show only the specified model. If missing, AIC is used
    to select the model. The model names (which can be abbreviated)
    are \code{Null}, \code{Preemption}, \code{Lognormal}, \code{Zipf},
    \code{Mandelbrot}. }
  \item{legend}{Add legend of line colours.}
  \item{as.table}{Arrange panels starting from upper left corner (passed
    to \code{\link[lattice]{xyplot}}).}
  \item{xlab,ylab}{Labels for \code{x} and \code{y} axes.}
  \item{type}{Type of the plot, \code{"b"} for plotting both observed points
    and fitted lines, \code{"p"} for only points, \code{"l"} for only
    fitted lines, and \code{"n"} for only setting the frame. }
  \item{log}{Use logarithmic scale for given axis. The default
    \code{log = "y"} gives the traditional plot of community ecology
    where the pre-emption model is a straight line, and with
    \code{log = "xy"} Zipf model is a straight line. With
    \code{log = ""} both axes are in the original arithmetic scale.}
  \item{\dots}{Other parameters to functions. }
}
\details{

  Rank--Abundance Dominance (RAD) or Dominance/Diversity plots
  (Whittaker 1965) display logarithmic species abundances against
  species rank order. These plots are supposed to be effective in
  analysing types of abundance distributions in communities. These
  functions fit some of the most popular models mainly following
  Wilson (1991).

  Functions \code{rad.null}, \code{rad.preempt}, \code{rad.lognormal},
  \code{rad.zipf} and \code{zipfbrot} fit the individual models
  (described below) for a single vector (row of data frame), and
  function \code{radfit} fits all models. The argument of the function
  \code{radfit} can be either a vector for a single community or a data
  frame where each row represents a distinct community.
  
  Function \code{rad.null} fits a brokenstick model where the expected
  abundance of species at rank \eqn{r} is \eqn{a_r = (J/S)
  \sum_{x=r}^S (1/x)}{a[r] = J/S sum(from x=r to S) 1/x} (Pielou
  1975), where \eqn{J} is the total number of individuals (site total)
  and \eqn{S} is the total number of species in the community.  This
  gives a Null model where the individuals are randomly distributed
  among observed species, and there are no fitted parameters.
  Function \code{rad.preempt} fits the niche preemption model,
  a.k.a. geometric series or Motomura model, where the expected
  abundance \eqn{a} of species at rank \eqn{r} is \eqn{a_r = J \alpha
  (1 - \alpha)^{r-1}}{a[r] = J*alpha*(1-alpha)^(r-1)}. The only
  estimated parameter is the preemption coefficient \eqn{\alpha} which
  gives the decay rate of abundance per rank.  The niche preemption
  model is a straight line in a RAD plot.  Function
  \code{rad.lognormal} fits a log-Normal model which assumes that the
  logarithmic abundances are distributed Normally, or \eqn{a_r = \exp(
  \log \mu + \log \sigma N)}{a[r] = exp(log(mu) + log(sigma) * N)},
  where \eqn{N} is a Normal deviate.  Function \code{rad.zipf} fits
  the Zipf model \eqn{a_r = J p_1 r^\gamma}{a[r] = J*p1*r^gamma} where
  \eqn{p_1}{p1} is the fitted proportion of the most abundant species,
  and \eqn{\gamma} is a decay coefficient. The Zipf--Mandelbrot model
  (\code{rad.zipfbrot}) adds one parameter: \eqn{a_r = J c (r +
  \beta)^\gamma}{a[r] = J*c*(r+beta)^gamma} after which \eqn{p_1}{p1}
  of the Zipf model changes into a meaningless scaling constant
  \eqn{c}. 

  Log-Normal and Zipf models are generalized linear models
  (\code{\link{glm}}) with logarithmic link function.  Zipf--Mandelbrot
  adds one nonlinear parameter to the Zipf model, and is fitted using
  \code{\link{nlm}} for the nonlinear parameter and estimating other
  parameters and log-Likelihood with \code{\link{glm}}. Preemption
  model is fitted as a purely nonlinear model. There are no estimated
  parameters in the Null model.  

  The default \code{\link{family}} is \code{poisson} which is
  appropriate only for genuine counts (integers), but other families
  that accept \code{link = "log"} can be used. Families
  \code{\link{Gamma}} or \code{\link{gaussian}} may be appropriate for
  abundance data, such as cover. The ``best'' model is selected by
  \code{\link{AIC}}. Therefore ``quasi'' families such as
  \code{\link{quasipoisson}} cannot be used: they do not have
  \code{\link{AIC}} nor log-Likelihood needed in non-linear models.
  
  All these functions have their own \code{plot} functions. When
  \code{radfit} was applied for a data frame, \code{plot} uses
  \code{\link[lattice]{Lattice}} graphics, and other \code{plot}
  functions use ordinary graphics. The ordinary graphics functions
  return invisibly an \code{\link{ordiplot}} object for observed points,
  and function \code{\link{identify.ordiplot}} can be used to label
  selected species.  Alternatively, \code{radlattice} uses
  \code{\link[lattice]{Lattice}} graphics to display each \code{radfit}
  model of a single site in a separate panel together with their AIC or
  BIC values.

  Function \code{as.rad} is a base function to construct ordered RAD
  data. Its \code{plot} is used by other RAD \code{plot} functions
  which pass extra arguments (such as \code{xlab} and \code{log}) to
  this function.

}

\value{

  Functions \code{rad.null}, \code{rad.preempt}, \code{rad.lognormal},
  \code{zipf} and \code{zipfbrot} fit each a single RAD model to a
  single site. The result object has class \code{"radline"} and
  inherits from \code{\link{glm}}, and can be handled by some (but not
  all) \code{\link{glm}} methods.

  Function \code{radfit} fits all models either to a single site or to
  all rows of a data frame or a matrix. When fitted to a single site,
  the function returns an object of class \code{"radfit"} with items
  \code{y} (observed values), \code{\link{family}}, and \code{models}
  which is a list of fitted \code{"radline"} models.  When applied for a
  data frame or matrix, \code{radfit} function returns an object of
  class \code{"radfit.frame"} which is a list of \code{"radfit"}
  objects, each item names by the corresponding row name.

  All result objects (\code{"radline"}, \code{"radfit"},
  \code{"radfit.frame"}) can be accessed with same method functions.
  The following methods are available: \code{\link{AIC}},
  \code{\link{coef}}, \code{\link{deviance}}, \code{\link{logLik}}. In
  addition the fit results can be accessed with \code{\link{fitted}},
  \code{\link{predict}} and \code{\link{residuals}} (inheriting from
  \code{\link{residuals.glm}}).The graphical functions were discussed
  above in Details.

}

\references{
  Pielou, E.C. (1975) \emph{Ecological Diversity}. Wiley & Sons.
  
  Preston, F.W. (1948) The commonness and rarity of
  species. \emph{Ecology} 29, 254--283.
  
  Whittaker, R. H. (1965) Dominance and diversity in plant
  communities. \emph{Science} 147, 250--260.

  Wilson, J. B. (1991) Methods for fitting dominance/diversity
  curves. \emph{Journal of Vegetation Science} 2, 35--46.
}
\author{ Jari Oksanen }
\note{
  The RAD models are usually fitted for proportions instead of original
  abundances. However, nothing in these models seems to require division
  of abundances by site totals, and original observations are used in
  these functions. If you wish to use proportions, you must standardize
  your data by site totals, e.g. with \code{\link{decostand}} and use
  appropriate \code{\link{family}} such as \code{\link{Gamma}}.

  The lognormal model is fitted in a standard way, but I do think this is
  not quite correct -- at least it is not equivalent to fitting Normal
  density to log abundances like originally suggested (Preston 1948). 
  
  Some models may fail. In particular, estimation of the Zipf-Mandelbrot
  model is difficult.  If the fitting fails, \code{NA} is returned.

  Wilson (1991) defined preemption model as \eqn{a_r = J p_1 (1
    - \alpha)^{r-1}}{a[r] = J*p1*(1 - alpha)^(r-1)}, where \eqn{p_1}{p1}
    is the fitted proportion of the first species. However, parameter
    \eqn{p_1}{p1} is completely defined by \eqn{\alpha} since the fitted
    proportions must add to one, and therefore I handle preemption as a
    one-parameter model.

    Veiled log-Normal model was included in earlier releases of this
    function, but it was removed because it was flawed: an implicit veil
    line also appears in the ordinary log-Normal. The latest release version
    with \code{rad.veil} was \code{1.6-10}. 
}
\seealso{\code{\link{fisherfit}} and \code{\link{prestonfit}}.
  An alternative approach is to use
  \code{\link{qqnorm}} or  \code{\link{qqplot}} with any distribution.
  For controlling graphics: \code{\link[lattice]{Lattice}},
  \code{\link[lattice]{xyplot}}, \code{\link[lattice]{lset}}. }
\examples{
data(BCI)
mod <- rad.lognormal(BCI[5,])
mod
plot(mod)
mod <- radfit(BCI[1,])
## Standard plot overlaid for all models
## Pre-emption model is a line
plot(mod)
## log for both axes: Zipf model is a line
plot(mod, log = "xy")
## Lattice graphics separately for each model
radlattice(mod)
# Take a subset of BCI to save time and nerves
mod <- radfit(BCI[3:5,])
mod
plot(mod, pch=".")
}
\keyword{ univar }
\keyword{ distribution }

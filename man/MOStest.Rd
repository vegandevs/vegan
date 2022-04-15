\encoding{UTF-8}
\name{MOStest}
\alias{MOStest}
\alias{plot.MOStest}
\alias{fieller.MOStest}
\alias{profile.MOStest}
\alias{confint.MOStest}

\title{ Mitchell-Olds \& Shaw Test for the Location of Quadratic Extreme }
\description{
  Mitchell-Olds & Shaw test concerns the location of the highest (hump)
  or lowest (pit) value of a quadratic curve at given points. Typically,
  it is used to study whether the quadratic hump or pit is located
  within a studied interval. The current test is generalized so that it
  applies generalized linear models (\code{\link{glm}}) with link
  function instead of simple quadratic curve.  The test was popularized
  in ecology for the analysis of humped species richness patterns
  (Mittelbach et al. 2001), but it is more general. With logarithmic
  link function, the quadratic response defines the Gaussian response
  model of ecological gradients (ter Braak & Looman 1986), and the test
  can be used for inspecting the location of Gaussian optimum within a
  given range of the gradient. It can also be used to replace Tokeshi's
  test of \dQuote{bimodal} species frequency distribution. 
}
\usage{
MOStest(x, y, interval, ...)
\method{plot}{MOStest}(x, which = c(1,2,3,6), ...)
fieller.MOStest(object, level = 0.95)
\method{profile}{MOStest}(fitted, alpha = 0.01, maxsteps = 10, del = zmax/5, ...)
\method{confint}{MOStest}(object, parm = 1, level = 0.95, ...)
}

\arguments{
  \item{x}{The independent variable or plotting object in \code{plot}. }
  \item{y}{The dependent variable. }
  \item{interval}{The two points at which the test statistic is
    evaluated. If missing, the extremes of \code{x} are used. }
  \item{which}{Subset of plots produced. Values \code{which=1} and
    \code{2} define plots specific to \code{MOStest} (see Details), and
    larger values select graphs of \code{\link{plot.lm}} (minus 2). }
  \item{object, fitted}{A result object from \code{MOStest}.}
  \item{level}{The confidence level required.}
  \item{alpha}{Maximum significance level allowed.}
  \item{maxsteps}{Maximum number of steps in the profile.}
  \item{del}{A step length parameter for the profile (see code).}
  \item{parm}{Ignored.}
  \item{\dots}{ Other variables passed to functions. Function
    \code{MOStest} passes these to \code{\link{glm}} so that
    these can include \code{\link{family}}. The other functions pass
    these to underlying graphical functions. }
}

\details{

  The function fits a quadratic curve \eqn{\mu = b_0 + b_1 x + b_2
  x^2} with given \code{\link{family}} and link function.  If \eqn{b_2
  < 0}, this defines a unimodal curve with highest point at \eqn{u =
  -b_1/(2 b_2)} (ter Braak & Looman 1986). If \eqn{b_2 > 0}, the
  parabola has a minimum at \eqn{u} and the response is sometimes
  called \dQuote{bimodal}.  The null hypothesis is that the extreme
  point \eqn{u} is located within the interval given by points
  \eqn{p_1} and \eqn{p_2}. If the extreme point \eqn{u} is exactly at
  \eqn{p_1}, then \eqn{b_1 = 0} on shifted axis \eqn{x - p_1}.  In the
  test, origin of \code{x} is shifted to the values \eqn{p_1} and
  \eqn{p_2}, and the test statistic is based on the differences of
  deviances between the original model and model where the origin is
  forced to the given location using the standard
  \code{\link{anova.glm}} function (Oksanen et al. 2001).
  Mitchell-Olds & Shaw (1987) used the first degree coefficient with
  its significance as estimated by the \code{\link{summary.glm}}
  function.  This give identical results with Normal error, but for
  other error distributions it is preferable to use the test based on
  differences in deviances in fitted models.

  The test is often presented as a general test for the location of the
  hump, but it really is dependent on the quadratic fitted curve. If the
  hump is of different form than quadratic, the test may be
  insignificant.

  Because of strong assumptions in the test, you should use the support
  functions to inspect the fit. Function \code{plot(..., which=1)}
  displays the data points, fitted quadratic model, and its approximate
  95\% confidence intervals (2 times SE). Function \code{plot} with
  \code{which = 2} displays the approximate confidence interval of
  the polynomial coefficients, together with two lines indicating the
  combinations of the coefficients that produce the evaluated points of
  \code{x}. Moreover, the cross-hair shows the approximate confidence
  intervals for the polynomial coefficients ignoring their
  correlations. Higher values of \code{which} produce corresponding
  graphs from \code{\link{plot.lm}}. That is, you must add 2 to the
  value of \code{which} in \code{\link{plot.lm}}.

  Function \code{fieller.MOStest} approximates the confidence limits
  of the location of the extreme point (hump or pit) using Fieller's
  theorem following ter Braak & Looman (1986). The test is based on
  quasideviance except if the \code{\link{family}} is \code{poisson}
  or \code{binomial}. Function \code{profile} evaluates the profile
  deviance of the fitted model, and \code{confint} finds the profile
  based confidence limits following Oksanen et al. (2001).

  The test is typically used in assessing the significance of diversity
  hump against productivity gradient (Mittelbach et al. 2001). It also
  can be used for the location of the pit (deepest points) instead of
  the Tokeshi test. Further, it can be used to test the location of the
  the Gaussian optimum in ecological gradient analysis (ter Braak &
  Looman 1986, Oksanen et al. 2001).
}

\value{
  The function is based on \code{\link{glm}}, and it returns the result
  of object of \code{glm} amended with the result of the test. The new
  items in the \code{MOStest} are: 
  \item{isHump }{\code{TRUE} if the response is a
    hump.}
  \item{isBracketed}{\code{TRUE} if the hump or the pit is bracketed by
    the evaluated points.} 
  \item{hump}{Sorted vector of location of the hump or the pit and the
    points where the test was evaluated.}
  \item{coefficients}{Table of test statistics and their significances.}
}

\references{
Mitchell-Olds, T. & Shaw, R.G. 1987. Regression analysis of natural
selection: statistical inference and biological
interpretation. \emph{Evolution} 41, 1149--1161.

Mittelbach, G.C. Steiner, C.F., Scheiner, S.M., Gross, K.L., Reynolds,
H.L., Waide, R.B., Willig, R.M., Dodson, S.I. & Gough, L. 2001. What is
the observed relationship between species richness and productivity?
\emph{Ecology} 82, 2381--2396.

Oksanen, J., Läärä, E., Tolonen, K. & Warner, B.G. 2001. Confidence
intervals for the optimum in the Gaussian response
function. \emph{Ecology} 82, 1191--1197.

ter Braak, C.J.F & Looman, C.W.N 1986. Weighted averaging, logistic
regression and the Gaussian response model. \emph{Vegetatio} 65,
3--11. 
}
\author{Jari Oksanen }

\note{ 
Function \code{fieller.MOStest} is based on package \pkg{optgrad} in
the Ecological Archives
(\url{https://figshare.com/articles/dataset/Full_Archive/3521975})
accompanying Oksanen et al. (2001). The Ecological Archive package
\pkg{optgrad} also contains profile deviance method for the location
of the hump or pit, but the current implementation of \code{profile}
and \code{confint} rather follow the example of
\code{\link[MASS]{profile.glm}} and \code{\link[MASS]{confint.glm}} in
the \pkg{MASS} package.
}

\seealso{The no-interaction model can be fitted with \code{\link{humpfit}}. }
\examples{
## The Al-Mufti data analysed in humpfit():
mass <- c(140,230,310,310,400,510,610,670,860,900,1050,1160,1900,2480)
spno <- c(1,  4,  3,  9, 18, 30, 20, 14,  3,  2,  3,  2,  5,  2)
mod <- MOStest(mass, spno)
## Insignificant
mod
## ... but inadequate shape of the curve
op <- par(mfrow=c(2,2), mar=c(4,4,1,1)+.1)
plot(mod)
## Looks rather like log-link with Poisson error and logarithmic biomass
mod <- MOStest(log(mass), spno, family=quasipoisson)
mod
plot(mod)
par(op)
## Confidence Limits
fieller.MOStest(mod)
confint(mod)
plot(profile(mod))
}

\keyword{ models }
\keyword{ regression }

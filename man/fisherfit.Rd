\name{fisherfit}
\alias{fisherfit}
\alias{as.fisher}
\alias{plot.fisherfit}
\alias{profile.fisherfit}
\alias{confint.fisherfit}
\alias{plot.profile.fisherfit}
\alias{prestonfit}
\alias{prestondistr}
\alias{as.preston}
\alias{plot.prestonfit}
\alias{lines.prestonfit}
\alias{plot.preston}
\alias{lines.preston}
\alias{plot.fisher}
\alias{veiledspec}

\title{Fit Fisher's Logseries and Preston's Lognormal Model to Abundance Data}
\description{
  Function \code{fisherfit} fits Fisher's logseries to abundance
  data. Function \code{prestonfit} groups species frequencies into
  doubling octave classes and fits Preston's lognormal model, and
  function \code{prestondistr} fits the truncated lognormal model
  without pooling the data into octaves.
}
\usage{
fisherfit(x, ...)
\method{confint}{fisherfit}(object, parm, level = 0.95, ...)
\method{profile}{fisherfit}(fitted, alpha = 0.01, maxsteps = 20, del = zmax/5, 
    ...)
prestonfit(x, tiesplit = TRUE, ...)
prestondistr(x, truncate = -1, ...)
\method{plot}{prestonfit}(x, xlab = "Frequency", ylab = "Species", bar.col = "skyblue", 
    line.col = "red", lwd = 2, ...)
\method{lines}{prestonfit}(x, line.col = "red", lwd = 2, ...)
veiledspec(x, ...)
as.fisher(x, ...)
\method{plot}{fisher}(x, xlab = "Frequency", ylab = "Species", bar.col = "skyblue",
             kind = c("bar", "hiplot", "points", "lines"), add = FALSE, ...)
as.preston(x, tiesplit = TRUE, ...)
\method{plot}{preston}(x, xlab = "Frequency", ylab = "Species", bar.col = "skyblue", ...)
\method{lines}{preston}(x, xadjust = 0.5, ...)
}

\arguments{
  \item{x}{Community data vector for fitting functions or their result
    object for \code{plot} functions.}
  \item{object, fitted}{Fitted model.}
  \item{parm}{Not used.}
  \item{level}{The confidence level required.}
  \item{alpha}{The extend of profiling as significance.}
  \item{maxsteps}{Maximum number of steps in profiling.}
  \item{del}{Step length.}
  \item{tiesplit}{Split frequencies \eqn{1, 2, 4, 8} etc between adjacent 
    octaves.}
  \item{truncate}{Truncation point for log-Normal model, in log2
    units. Default value \eqn{-1} corresponds to the left border of zero
    Octave. The choice strongly influences the fitting results.}
  \item{xlab, ylab}{Labels for \code{x} and \code{y} axes.}
  \item{bar.col}{Colour of data bars.}
  \item{line.col}{Colour of fitted line.}
  \item{lwd}{Width of fitted line.}
  \item{kind}{Kind of plot to drawn: \code{"bar"} is similar bar plot
    as in \code{plot.fisherfit}, \code{"hiplot"} draws vertical lines
    as with \command{plot(..., type="h")}, and \code{"points"} and
    \code{"lines"} are obvious.}
  \item{add}{Add to an existing plot.}
  \item{xadjust}{Adjustment of horizontal positions in octaves.}
  \item{\dots}{Other parameters passed to functions. Ignored in 
    \code{prestonfit} and \code{tiesplit} passed to \code{as.preston} in
    \code{prestondistr}. }
}
\details{
  In Fisher's logarithmic series the expected
  number of species \eqn{f} with \eqn{n} observed individuals is
  \eqn{f_n = \alpha x^n / n} (Fisher et al. 1943). The estimation
  follows Kempton & Taylor (1974) and uses function
  \code{\link{nlm}}. The estimation is possible only for genuine
  counts of individuals. The parameter \eqn{\alpha} is used as a
  diversity index, and \eqn{\alpha} and its standard error can be
  estimated with a separate function \code{\link{fisher.alpha}}. The
  parameter \eqn{x} is taken as a nuisance parameter which is not
  estimated separately but taken to be \eqn{n/(n+\alpha)}. Helper
  function \code{as.fisher} transforms abundance data into Fisher
  frequency table.

  Function  \code{fisherfit} estimates the standard error of
  \eqn{\alpha}{alpha}. However, the confidence limits cannot be directly
  estimated from the standard errors, but you should use function
  \code{confint} based on profile likelihood. Function \code{confint}
  uses function \code{\link[MASS]{confint.glm}} of the \pkg{MASS}
  package, using \code{profile.fisherfit} for the profile
  likelihood. Function \code{profile.fisherfit} follows
  \code{\link[MASS]{profile.glm}} and finds the \eqn{\tau}{tau} parameter or
  signed square root of two times log-Likelihood profile. The profile can
  be inspected with a \code{plot} function which shows the \eqn{\tau}{tau}
  and a dotted line corresponding to the Normal assumption: if standard
  errors can be directly used in Normal inference these two lines
  are similar.

  Preston (1948) was not satisfied with Fisher's model which seemed to
  imply infinite species richness, and postulated that rare species is
  a diminishing class and most species are in the middle of frequency
  scale. This was achieved by collapsing higher frequency classes into
  wider and wider ``octaves'' of doubling class limits: 1, 2, 3--4,
  5--8, 9--16 etc. occurrences. It seems that Preston regarded
  frequencies 1, 2, 4, \emph{etc.}. as \dQuote{tied} between octaves
  (Williamson & Gaston 2005). This means that only half of the species
  with frequency 1 are shown in the lowest octave, and the rest are
  transferred to the second octave. Half of the species from the
  second octave are transferred to the higher one as well, but this is
  usually not as large a number of species. This practise makes data
  look more lognormal by reducing the usually high lowest
  octaves. This can be achieved by setting argument \code{tiesplit =
  TRUE}. With \code{tiesplit = FALSE} the frequencies are not split,
  but all ones are in the lowest octave, all twos in the second, etc.
  Williamson & Gaston (2005) discuss alternative definitions in
  detail, and they should be consulted for a critical review of
  log-Normal model.

  Any logseries data will look like lognormal when plotted in
  Preston's way. The expected frequency \eqn{f} at abundance octave
  \eqn{o} is defined by \eqn{f_o = S_0 \exp(-(\log_2(o) -
  \mu)^2/2/\sigma^2)}{f = S0 exp(-(log2(o)-mu)^2/2/sigma^2)}, where
  \eqn{\mu} is the location of the mode and \eqn{\sigma} the width,
  both in \eqn{\log_2}{log2} scale, and \eqn{S_0}{S0} is the expected
  number of species at mode. The lognormal model is usually truncated
  on the left so that some rare species are not observed. Function
  \code{prestonfit} fits the truncated lognormal model as a second
  degree log-polynomial to the octave pooled data using Poisson (when
  \eqn{tiesplit = FALSE}) or quasi-Poisson (when \eqn{tiesplit =
  TRUE}).  error. Function \code{prestondistr} fits left-truncated
  Normal distribution to \eqn{\log_2}{log2} transformed non-pooled
  observations with direct maximization of log-likelihood. Function
  \code{prestondistr} is modelled after function
  \code{\link[MASS]{fitdistr}} which can be used for alternative
  distribution models. 

  The functions have common \code{print}, \code{plot} and \code{lines}
  methods. The \code{lines} function adds the fitted curve to the
  octave range with line segments showing the location of the mode and
  the width (sd) of the response. Function \code{as.preston}
  transforms abundance data to octaves.  Argument \code{tiesplit} will
  not influence the fit in \code{prestondistr}, but it will influence
  the barplot of the octaves. 

  The total extrapolated richness from a fitted Preston model can be
  found with function \code{veiledspec}. The function accepts results
  both from \code{prestonfit} and from \code{prestondistr}. If
  \code{veiledspec} is called with a species count vector, it will
  internally use \code{prestonfit}. Function \code{\link{specpool}}
  provides alternative ways of estimating the number of unseen
  species. In fact, Preston's lognormal model seems to be truncated at
  both ends, and this may be the main reason why its result differ
  from lognormal models fitted in Rank--Abundance diagrams with
  functions \code{\link{rad.lognormal}}.  
}

\value{ The function \code{prestonfit} returns an object with fitted
  \code{coefficients}, and with observed (\code{freq}) and fitted
  (\code{fitted}) frequencies, and a string describing the fitting
  \code{method}. Function \code{prestondistr} omits the entry
  \code{fitted}.  The function \code{fisherfit} returns the result of
  \code{\link{nlm}}, where item \code{estimate} is \eqn{\alpha}. The
  result object is amended with the following items:
  \item{df.residuals}{Residual degrees of freedom.}
  \item{nuisance}{Parameter \eqn{x}.}  \item{fisher}{Observed data
  from \code{as.fisher}.}

}
\references{
  Fisher, R.A., Corbet, A.S. & Williams, C.B. (1943). The relation
  between the number of species and the number of individuals in a
  random sample of animal population. \emph{Journal of Animal Ecology}
  12: 42--58.

  Kempton, R.A. & Taylor, L.R. (1974). Log-series and log-normal
  parameters as diversity discriminators for
  Lepidoptera. \emph{Journal of Animal Ecology} 43: 381--399.

  Preston, F.W. (1948) The commonness and rarity of
  species. \emph{Ecology} 29, 254--283.

  Williamson, M. & Gaston, K.J. (2005). The lognormal distribution is
  not an appropriate null hypothesis for the species--abundance
  distribution. \emph{Journal of Animal Ecology} 74, 409--422.
}

\author{Bob O'Hara \email{bob.ohara@helsinki.fi} (\code{fisherfit})
and Jari Oksanen. }

\seealso{\code{\link{diversity}}, \code{\link{fisher.alpha}},
  \code{\link{radfit}}, \code{\link{specpool}}. Function
  \code{\link[MASS]{fitdistr}} of \pkg{MASS} package was used as the
  model for \code{prestondistr}. Function \code{\link{density}} can be used for
  smoothed ``non-parametric'' estimation of responses, and
  \code{\link{qqplot}} is an alternative, traditional and more effective
  way of studying concordance of observed abundances to any distribution model.
}
\examples{
data(BCI)
mod <- fisherfit(BCI[5,])
mod
plot(profile(mod))
confint(mod)
# prestonfit seems to need large samples
mod.oct <- prestonfit(colSums(BCI))
mod.ll <- prestondistr(colSums(BCI))
mod.oct
mod.ll
plot(mod.oct)  
lines(mod.ll, line.col="blue3") # Different
## Smoothed density
den <- density(log2(colSums(BCI)))
lines(den$x, ncol(BCI)*den$y, lwd=2) # Fairly similar to mod.oct
## Extrapolated richness
veiledspec(mod.oct)
veiledspec(mod.ll)
}
\keyword{ univar }
\keyword{ distribution }

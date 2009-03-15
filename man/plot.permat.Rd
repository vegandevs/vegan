\encoding{UTF-8}
\name{diagnose.permat}
\alias{plot.permat}
\alias{lines.permat}
\alias{as.ts.permat}
\alias{as.mcmc.permat}

\title{Plotting Methods for Matrix Permutation Algorithms}

\description{ Plotting and time series extraction methods for 'permat' objects, that are useful in evaluating the performance of different null model settings, especially that of the sequential algorithms.}

\usage{
\method{plot}{permat}(x, type = "bray", ylab, xlab, col, lty,
    lowess = TRUE, plot = TRUE, text = TRUE, ...)
\method{lines}{permat}(x, type = "bray", ...)
\method{as.ts}{permat}(x, type = "bray", ...)
\method{as.mcmc}{permat}(x)
}
\arguments{
  \item{x}{Object of class \code{"permat"} created by the functions \code{\link{permatfull}} and \code{\link{permatswap}}.}
  \item{ylab, xlab, col, lty}{graphical parameters for the \code{plot} method.}
  \item{type}{Character, type of plot to be displayed: \code{"bray"} for Bray-Curtis dissimilarities, \code{"chisq"} for Chi-squared values.}
  \item{lowess, plot, text}{Logical arguments for the \code{plot} method, whether a lowess curve should be drawn, the plot should be drawn, and statistic values should be printed on the plot.}
  \item{\dots}{Other arguments passed to methods.}
} 

\details{

The \code{plot} method is useful for graphically testing for trend and independence of permuted matrices. This is especially important when using sequential algorithms (\code{"swap", "tswap", "abuswap"}).


The \code{as.ts} method can be used to extract Bray-Curtis dissimilarities or Chi-squared values as time series. This can further used in testing independence (see Examples). The method \code{as.mcmc} is useful for accessing diagnostic tool available in the 'code' package.
}

\value{
The \code{plot} creates aplor as a side effect.

The \code{as.ts} method returns an object of class 'ts'.
}
\author{\enc{P\'eter S\'olymos}{Peter Solymos}, \email{solymos@ualberta.ca}}

\seealso{
For functions to create 'permat' objects: \code{\link{permatfull}}, \code{\link{permatswap}}

For time-series diagnostics: \code{\link{Box.test}}, \code{\link{lag.plot}}, \code{\link{tsdiag}}, \code{\link{ar}}, \code{\link{arima}}
}
\examples{
data(BCI)
## Not sequential algorithm
a <- permatswap(BCI, "quasiswap")
## Sequential algorithm
b <- permatswap(BCI, "abuswap", fixedmar="col",
burnin=0, thin=100, times=50)
opar <- par(mfrow=c(2,2))
plot(a, main="Not sequential")
plot(b, main="Sequential")
plot(a, "chisq")
plot(b, "chisq")
par(opar)
## Extract Bray-Curtis dissimilarities
## as time series
bc <- as.ts(b)
## Lag plot
lag.plot(bc)
## First order autoregressive model
mar <- arima(bc, c(1,0,0))
mar
## Ljung-Box test of residuals
Box.test(mar$residuals)
## Graphical diagnostics
tsdiag(mar)
}

\keyword{multivariate}

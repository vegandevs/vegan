\encoding{UTF-8}
\name{diagnose.permat}
\alias{plot.permat}
\alias{lines.permat}
\alias{lagplot.permat}
\alias{lagplot}
\alias{diagnose.permat}
\alias{print.diagnose.permat}
\alias{diagnose}

\title{Diagnostic Tool for Matrix Permutation Algorithms}

\description{ Plotting and diagnostic methods for 'permat' objects, that are useful in evaluating the performance of different null model settings, especially that of the sequential algorithms.}

\usage{
\method{plot}{permat}(x, type = "bray", ylab, xlab, col, lty,
    lowess = TRUE, plot = TRUE, text = TRUE, ...)
\method{lines}{permat}(x, type = "bray", ...)
lagplot(x, ...)
\method{lagplot}{permat}(x, type = "bray", ...)
diagnose(x, ...)
\method{diagnose}{permat}(x, type = "bray", ...)
\method{print}{diagnose.permat}(x, ...)
}
\arguments{
  \item{x}{Object of class \code{"permat"} created by the functions \code{\link{permatfull}} and \code{\link{permatswap}}.}
  \item{ylab, xlab, col, lty}{graphical parameters for the \code{plot} method.}
  \item{type}{Character, type of plot to be displayed: \code{"bray"} for Bray-Curtis dissimilarities, \code{"chisq"} for Chi-squared values.}
  \item{lowess, plot, text}{Logical arguments for the \code{plot} method, whether a lowess curve should be drawn, the plot should be drawn, and statistic values should be printed on the plot.}
  \item{\dots}{Other arguments passed to methods.}
} 

\details{

The functions \code{lagplot} \code{diagnose} methods are useful for testing independence of permuted matrices
when using sequential algorithms (\code{"swap", "tswap", "abuswap"}).

}

\value{
The \code{diagnose} returns a list including first order autoregressive model and Ljung-Box test results to diagnose independence of the statistics (Bray-Curtis, Chi-square) derived from successive permuted matrices.
}
\author{\enc{P\'eter S\'olymos}{Peter Solymos}, \email{solymos@ualberta.ca} and Jari Oksanen}

\seealso{
For functions to create 'permat' objects: \code{\link{permatfull}}, \code{\link{permatswap}}

For time-series diagnostics: \code{\link{Box.test}}, \code{\link{lag.plot}}, \code{\link{tsdiag}}
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
## The lagplot is useful to diagnose independence
## of successive simulations for sequential algorithms
lagplot.permat(b) # highly non-independent
## Time series diagnostics
diagnose(b)
}

\keyword{multivariate}

\encoding{UTF-8}
\name{permat}
\alias{permatfull}
\alias{permatswap}
\alias{summary.permat}
\alias{print.summary.permat}
\alias{print.permat}
\alias{plot.permat}
\alias{lines.permat}
\alias{as.ts.permat}
\alias{as.mcmc.permat}

\title{Matrix Permutation Algorithms for Presence-Absence and Count Data}

\description{ Individual (for count data) or incidence (for
presence-absence data) based null models can be generated for
community level simulations. Options for preserving characteristics of
the original matrix (rows/columns sums, matrix fill) and
restricted permutations (based on strata) are discussed in the
Details section.}

\usage{
permatfull(m, fixedmar = "both", shuffle = "both", strata = NULL, 
    mtype = "count", times = 99, ...)
permatswap(m, method = "quasiswap", fixedmar="both", shuffle = "both",
    strata = NULL, mtype = "count", times = 99, 
    burnin = 0, thin = 1, ...)
\method{print}{permat}(x, digits = 3, ...)
\method{summary}{permat}(object, ...)
\method{print}{summary.permat}(x, digits = 2, ...)
\method{plot}{permat}(x, type = "bray", ylab, xlab, col, lty,
    lowess = TRUE, plot = TRUE, text = TRUE, ...)
\method{lines}{permat}(x, type = "bray", ...)
\method{as.ts}{permat}(x, type = "bray", ...)
\method{as.mcmc}{permat}(x)
}
\arguments{
  \item{m}{A community data matrix with plots (samples) as rows and
    species (taxa) as columns.} 
  \item{fixedmar}{character, stating which of the row/column sums should
    be preserved (\code{"none", "rows", "columns", "both"}).} 
  \item{strata}{Numeric vector or factor with length same as
    \code{nrow(m)} for grouping rows within strata for restricted
    permutations. Unique values or levels are used.} 
  \item{mtype}{Matrix data type, either \code{"count"} for count data,
    or \code{"prab"} for presence-absence type incidence data.} 
  \item{times}{Number of permuted matrices.} 
  \item{method}{Character for method used for the swap algorithm
    (\code{"swap"}, \code{"tswap"}, \code{"quasiswap"},
    \code{"backtrack"}) as described for function
    \code{\link{make.commsim}}. If \code{mtype="count"} the
    \code{"quasiswap"}, \code{"swap"}, \code{"swsh"} and
    \code{"abuswap"} methods are available (see details).} 
  \item{shuffle}{Character, indicating whether individuals
    (\code{"ind"}), samples (\code{"samp"}) or both (\code{"both"})
    should be shuffled, see details.} 
  \item{burnin}{Number of null communities discarded before proper
    analysis in sequential (\code{"swap", "tswap"}) methods.} 
  \item{thin}{Number of discarded permuted matrices between two
    evaluations in sequential (\code{"swap", "tswap"}) methods.} 
  \item{x, object}{Object of class \code{"permat"}} 
  \item{digits}{Number of digits used for rounding.}
  \item{ylab, xlab, col, lty}{graphical parameters for the \code{plot}
    method.} 
  \item{type}{Character, type of plot to be displayed: \code{"bray"} for
    Bray-Curtis dissimilarities, \code{"chisq"} for Chi-squared values.} 
  \item{lowess, plot, text}{Logical arguments for the \code{plot}
    method, whether a locally weighted regression curve should be drawn,
    the plot should be drawn, and statistic values should be printed on
    the plot.} 
  \item{\dots}{Other arguments passed to \code{\link{simulate.nullmodel}} 
    or methods.}
}

\details{
  The function \code{permatfull} is useful when matrix fill is
  allowed to vary, and matrix type is \code{count}.  The \code{fixedmar}
  argument is used to set constraints for permutation.  If \code{none}
  of the margins are fixed, cells are randomised within the matrix.  If
  \code{rows} or \code{columns} are fixed, cells within rows or columns
  are randomised, respectively.  If \code{both} margins are fixed, the
  \code{\link{r2dtable}} function is used that is based on Patefield's
  (1981) algorithm. For presence absence data, matrix fill should be
  necessarily fixed, and \code{permatfull} is a wrapper for the function
  \code{\link{make.commsim}}. The \code{r00, r0, c0, quasiswap}
  algorithms of \code{\link{make.commsim}} are used for \code{"none",
  "rows", "columns", "both"} values of the \code{fixedmar} argument,
  respectively

  The \code{shuffle} argument only have effect if the \code{mtype =
  "count"} and \code{permatfull} function is used with \code{"none",
  "rows", "columns"} values of \code{fixedmar}. All other cases for
  count data are individual based randomisations. The \code{"samp"} and
  \code{"both"} options result fixed matrix fill. The \code{"both"}
  option means that individuals are shuffled among non zero cells
  ensuring that there are no cell with zeros as a result, then cell
  (zero and new valued cells) are shuffled.

  The function \code{permatswap} is useful when with matrix fill
  (i.e. the proportion of empty cells) and row/columns sums should be
  kept constant. \code{permatswap} uses different kinds of swap
  algorithms, and row and columns sums are fixed in all cases.  For
  presence-absence data, the \code{swap} and \code{tswap} methods of
  \code{\link{make.commsim}} can be used.  For count data, a special
  swap algorithm ('swapcount') is implemented that results in permuted
  matrices with fixed marginals and matrix fill at the same time.

  The 'quasiswapcount' algorithm (\code{method="quasiswap"} and
  \code{mtype="count"}) uses the same trick as Carsten Dormann's
  \code{swap.web} function in the package
  \pkg{bipartite}. First, a random matrix is generated by the
  \code{\link{r2dtable}} function retaining row and column sums. Then
  the original matrix fill is reconstructed by sequential steps to
  increase or decrease matrix fill in the random matrix. These steps are
  based on swapping 2x2 submatrices (see 'swapcount' algorithm for
  details) to maintain row and column totals. This algorithm generates
  independent matrices in each step, so \code{burnin} and \code{thin}
  arguments are not considered. This is the default method, because this
  is not sequential (as \code{swapcount} is) so independence of subsequent
  matrices does not have to be checked.

  The \code{swapcount} algorithm (\code{method="swap"} and
  \code{mtype="count"}) tries to find 2x2 submatrices (identified by 2
  random row and 2 random column indices), that can be swapped in order
  to leave column and row totals and fill unchanged. First, the
  algorithm finds the largest value in the submatrix that can be swapped
  (\eqn{d}) and whether in diagonal or antidiagonal way. Submatrices
  that contain values larger than zero in either diagonal or
  antidiagonal position can be swapped. Swap means that the values in
  diagonal or antidiagonal positions are decreased by \eqn{d}, while
  remaining cells are increased by \eqn{d}. A swap is made only if fill
  doesn't change. This algorithm is sequential, subsequent matrices are
  not independent, because swaps modify little if the matrix is
  large. In these cases many burnin steps and thinning is needed to get
  independent random matrices. Although this algorithm is implemented in
  C, large burnin and thin values can slow it down
  considerably. WARNING: according to simulations, this algorithm seems
  to be biased and non random, thus its use should be avoided!

  The algorithm \code{"swsh"} in the function \code{permatswap} is a
  hybrid algorithm. First, it makes binary quasiswaps to keep row and
  column incidences constant, then non-zero values are modified
  according to the \code{shuffle} argument (only \code{"samp"} and
  \code{"both"} are available in this case, because it is applied only
  on non-zero values). It also recognizes the \code{fixedmar}
  argument which cannot be \code{"both"} (\pkg{vegan} versions <= 2.0
  had this algorithm with \code{fixedmar = "none"}).

  The algorithm \code{"abuswap"} produces two kinds of null models
  (based on \code{fixedmar="columns"} or \code{fixedmar="rows"}) as
  described in Hardy (2008; randomization scheme 2x and 3x,
  respectively).  These preserve column and row occurrences, and column
  or row sums at the same time. (Note that similar constraints
  can be achieved by the non sequential \code{"swsh"} algorithm
  with \code{fixedmar} argument set to \code{"columns"} or
  \code{"rows"}, respectively.)

  Constraints on row/column sums, matrix fill, total sum and sums within
  strata can be checked by the \code{summary} method. \code{plot} method
  is for visually testing the randomness of the permuted matrices,
  especially for the sequential swap algorithms. If there are any
  tendency in the graph, higher \code{burnin} and \code{thin} values can
  help for sequential methods.  New lines can be added to existing plot
  with the \code{lines} method.

  Unrestricted and restricted permutations: if \code{strata} is
  \code{NULL}, functions perform unrestricted permutations. Otherwise,
  it is used for restricted permutations. Each strata should contain at
  least 2 rows in order to perform randomization (in case of low row
  numbers, swap algorithms can be rather slow). If the design is not
  well balanced (i.e. same number of observations within each stratum),
  permuted matrices may be biased because same constraints are forced on
  submatrices of different dimensions. This often means, that the number
  of potential permutations will decrease with their dimensions.  So the
  more constraints we put, the less randomness can be expected.

  The \code{plot} method is useful for graphically testing for trend and
  independence of permuted matrices. This is especially important when
  using sequential algorithms (\code{"swap", "tswap", "abuswap"}).

  The \code{as.ts} method can be used to extract Bray-Curtis
  dissimilarities or Chi-squared values as time series. This can further
  used in testing independence (see Examples). The method \code{as.mcmc}
  is useful for accessing diagnostic tools available in the \pkg{coda}
  package.  }

\value{Functions \code{permatfull} and \code{permatswap} return an
  object of class \code{"permat"} containing the the function call
  (\code{call}), the original data matrix used for permutations
  (\code{orig}) and a list of permuted matrices with length \code{times}
  (\code{perm}).

  The \code{summary} method returns various statistics as a list
  (including mean Bray-Curtis dissimilarities calculated pairwise among
  original and permuted matrices, Chi-square statistics, and check
  results of the constraints; see Examples). Note that when
  \code{strata} is used in the original call, summary calculation may
  take longer.

  The \code{plot} creates a plot as a side effect.

  The \code{as.ts} method returns an object of class \code{"ts"}.  }


\references{ Original references for presence-absence algorithms are
  given on help page of \code{\link{make.commsim}}.

  Hardy, O. J. (2008) Testing the spatial phylogenetic structure of
  local communities: statistical performances of different null models
  and test statistics on a locally neutral community. Journal of Ecology
  96, 914--926. 

  Patefield, W. M. (1981) Algorithm AS159. An efficient method of
  generating r x c tables with given row and column totals.  
  Applied Statistics 30, 91--97.
}

\author{\enc{Péter Sólymos}{Peter Solymos},
\email{solymos@ualberta.ca} and Jari Oksanen}

\seealso{ For other functions to permute matrices:
\code{\link{make.commsim}}, \code{\link{r2dtable}},
\code{\link{sample}}.

For the use of these permutation algorithms: \code{\link{oecosimu}},
\code{\link{adipart}}, \code{\link{hiersimu}}.

For time-series diagnostics: \code{\link{Box.test}},
\code{\link{lag.plot}}, \code{\link{tsdiag}}, \code{\link{ar}},
\code{\link{arima}} 

For underlying low level implementation:
\code{\link{commsim}} and \code{\link{nullmodel}}.}

\examples{
## A simple artificial community data matrix.
m <- matrix(c(
    1,3,2,0,3,1,
    0,2,1,0,2,1,
    0,0,1,2,0,3,
    0,0,0,1,4,3
    ), 4, 6, byrow=TRUE)
## Using the quasiswap algorithm to create a 
## list of permuted matrices, where
## row/columns sums and matrix fill are preserved:
x1 <- permatswap(m, "quasiswap")
summary(x1)
## Unrestricted permutation retaining
## row/columns sums but not matrix fill:
x2 <- permatfull(m)
summary(x2)
## Unrestricted permutation of presence-absence type
## not retaining row/columns sums:
x3 <- permatfull(m, "none", mtype="prab")
x3$orig  ## note: original matrix is binarized!
summary(x3)
## Restricted permutation,
## check sums within strata:
x4 <- permatfull(m, strata=c(1,1,2,2))
summary(x4)

## NOTE: 'times' argument usually needs to be >= 99
## here much lower value is used for demonstration

## Not sequential algorithm
data(BCI)
a <- permatswap(BCI, "quasiswap", times=19)
## Sequential algorithm
b <- permatswap(BCI, "abuswap", fixedmar="col",
    burnin=0, thin=100, times=19)
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
Box.test(residuals(mar))
## Graphical diagnostics
tsdiag(mar)
}

\keyword{multivariate}
\keyword{datagen}

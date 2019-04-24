\name{permutations}
\alias{permutations}

\title{Permutation tests in Vegan}
\description{
  From version 2.2-0, \pkg{vegan} has significantly improved access to
  restricted permutations which brings it into line with those offered
  by Canoco. The permutation designs are modelled after the permutation
  schemes of Canoco 3.1 (ter Braak, 1990).

  \pkg{vegan} currently provides for the following features within
  permutation tests:
  \enumerate{
    \item{Free permutation of \emph{DATA}, also known as randomisation,}
    \item{Free permutation of \emph{DATA} within the levels of a
      grouping variable,}
    \item{Restricted permutations for line transects or time series,}
    \item{Permutation of groups of samples whilst retaining the
      within-group ordering,}
    \item{Restricted permutations for spatial grids,}
    \item{Blocking, samples are never permuted \emph{between} blocks,
      and}
    \item{Split-plot designs, with permutation of whole plots, split
      plots, or both.}
  }
  Above, we use \emph{DATA} to mean either the observed data themselves
  or some function of the data, for example the residuals of an
  ordination model in the presence of covariables.

  These capabilities are provided by functions from the \pkg{permute}
  package. The user can request a particular type of permutation by
  supplying the \code{permutations} argument of a function with an
  object returned by \code{\link{how}}, which defines how samples should
  be permuted. Alternatively, the user can simply specify the required
  number of permutations and a simple randomisation procedure will be
  performed. Finally, the user can supply a matrix of permutations (with
  number of rows equal to the number of permutations and number of
  columns equal to the number of observations in the data) and
  \pkg{vegan} will use these permutations instead of generating new
  permutations.

  The majority of functions in \pkg{vegan} allow for the full range of
  possibilities outlined above. Exceptions include
  \code{\link{kendall.post}} and \code{\link{kendall.global}}.

  The Null hypothesis for the first two types of permutation test listed
  above assumes free exchangeability of \emph{DATA} (within the levels
  of the grouping variable, if specified). Dependence between
  observations, such as that which arises due to spatial or temporal
  autocorrelation, or more-complicated experimental designs, such as
  split-plot designs, violates this fundamental assumption of the test
  and requires more complex restricted permutation test designs. It is
  these designs that are available via the \pkg{permute} package and to
  which \pkg{vegan} provides access from version 2.2-0 onwards.

  Unless otherwise stated in the help pages for specific functions,
  permutation tests in \pkg{vegan} all follow the same format/structure:
  \enumerate{
    \item{An appropriate test statistic is chosen. Which statistic is
      chosen should be described on the help pages for individual
      functions.}
    \item{The value of the test statistic is evaluate for the observed
      data and analysis/model and recorded. Denote this value
      \eqn{x_0}{x[0]}.}
    \item{The \emph{DATA} are randomly permuted according to one of the
      above schemes, and the value of the test statistic for this
      permutation is evaluated and recorded.}
    \item{Step 3 is repeated a total of \eqn{n} times, where \eqn{n} is
      the number of permutations requested. Denote these values as
      \eqn{x_i}{x[i]}, where \eqn{i = 1, ..., n}{{i = 1, \ldots, n}.}}
    \item{Count the number of values of the test statistic,
      \eqn{x_i}{x[i]}, in the Null distribution that are as extreme as
      test statistic for the observed data \eqn{x_0}{x[0]}. Denote this
      count as \eqn{N}.

      We use the phrase \emph{as extreme} to include cases where a
      two-sided test is performed and large negative values of the test
      statistic should be considered.}
    \item{The permutation p-value is computed as
      \deqn{p = \frac{N + 1}{n + 1}}{(N + 1) / (n + 1)}}
  }

  The above description illustrates why the default number of
  permutations specified in \pkg{vegan} functions takes values of 199 or
  999 for example. Pretty \emph{p} values are achieved because the
  \eqn{+ 1} in the denominator results in division by 200 or 1000, for
  the 199 or 999 random permutations used in the test.

  The simple intuition behind the presence of \eqn{+ 1} in the numerator
  and denominator is that these represent the inclusion of the observed
  value of the statistic in the Null distribution (e.g. Manly 2006).
  Phipson & Smyth (2010) present a more compelling explanation for the
  inclusion of \eqn{+ 1} in the numerator and denominator of the
  \emph{p} value calculation.

  Fisher (1935) had in mind that a permutation test would involve
  enumeration of all possible permutations of the data yielding an exact
  test. However, doing this complete enumeration may not be feasible in
  practice owing to the potentially vast number of arrangements of the
  data, even in modestly-sized data sets with free permutation of
  samples. As a result we evaluate the \emph{p} value as the tail
  probability of the Null distribution of the test statistic directly
  from the random sample of possible permutations. Phipson & Smyth
  (2010) show that the naive calculation of the permutation \emph{p}
  value is

  \deqn{p = \frac{N}{n}}{p = N / n}

  which leads to an invalid test with incorrect type I error rate. They
  go on to show that by replacing the unknown tail probability (the
  \emph{p} value) of the Null distribution with the biased estimator

  \deqn{p = \frac{N + 1}{n + 1}}{p = (N + 1) / (n + 1)}

  that the positive bias induced is of just the right size to
  account for the  uncertainty in the estimation of the tail probability
  from the set of randomly sampled permutations to yield a test with the
  correct type I error rate.

  The estimator described above is correct for the situation where
  permutations of the data are samples randomly \emph{without}
  replacement. This is not strictly what happens in \pkg{vegan} because
  permutations are drawn pseudo-randomly independent of one
  another. Note that the actual chance of this happening is practice is
  small but the functions in \pkg{permute} do not guarantee to generate
  a unique set of permutations unless complete enumeration of
  permutations is requested. This is not feasible for all but the
  smallest of data sets or restrictive of permutation designs, but in
  such cases the chance of drawing a set of permutations with repeats is
  lessened as the sample size, and thence the size of set of all
  possible permutations, increases.

  Under the situation of sampling permutations with replacement then,
  the tail probability \eqn{p} calculated from the biased estimator
  described above is somewhat \strong{conservative}, being too large by
  an amount that depends on the number of possible values that the test
  statistic can take under permutation of the data (Phipson & Smyth,
  2010). This represents a slight loss of statistical power for the
  conservative \emph{p} value calculation used here. However, unless
  sample sizes are small and the the permutation design such that the
  set of values that the test statistic can take is also small, this
  loss of power is unlikely to be critical.

  The minimum achievable p-value is

  \deqn{p_{\mathrm{min}} = \frac{1}{n + 1}}{p[min] = 1 / (n + 1)}

  and hence depends on the number of permutations evaluated. However,
  one cannot simply increase the number of permutations (\eqn{n}) to
  achieve a potentially lower p-value unless the number of observations
  available permits such a number of permutations. This is unlikely to
  be a problem for all but the smallest data sets when free permutation
  (randomisation) is valid, but in restricted permutation designs with a
  low number of observations, there may not be as many unique
  permutations of the data as you might desire to reach the required
  level of significance.

  It is currently the responsibility of the user to determine the total
  number of possible permutations for their \emph{DATA}. The number of
  possible permutations allowed under the specified design can be
  calculated using \code{\link[permute]{numPerms}} from the
  \pkg{permute} package. Heuristics employed within the
  \code{\link[permute]{shuffleSet}} function used by \pkg{vegan} can be
  triggered to generate the entire set of permutations instead of a
  random set. The settings controlling the triggering of the complete
  enumeration step are contained within a permutation design created
  using \code{link[permute]{how}} and can be set by the user. See
  \code{\link[permute]{how}} for details.

  Limits on the total number of permutations of \emph{DATA} are more
  severe in temporally or spatially ordered data or experimental designs
  with low replication. For example, a time series of \eqn{n = 100}
  observations has just 100 possible permutations \strong{including} the
  observed ordering.

  In situations where only a low number of permutations is possible due
  to the nature of \emph{DATA} or the experimental design, enumeration
  of all permutations becomes important and achievable computationally.

  Above, we have provided only a brief overview of the capabilities of
  \pkg{vegan} and \pkg{permute}. To get the best out of the new
  functionality and for details on how to set up permutation designs
  using \code{\link[permute]{how}}, consult the vignette
  \emph{Restricted permutations; using the permute package} supplied
  with \pkg{permute} and accessible via \code{vignette("permutations",
  package = "permute").}
}

\section{Random Number Generation}{

  The permutations are based on the random number generator provided
  by \R. This may change in \R{} releases and change the permutations
  and \pkg{vegan} test results. One such change was in \R{} release
  3.6.0. The new version is clearly better for permutation tests and
  you should use it. However, if you need to reproduce old results,
  you can set the \R{} random number generator to a previous version
  with \code{\link{RNGversion}}.

} % section RNG

\seealso{
  \code{\link{permutest}} for the main interface in \pkg{vegan}. See
  also \code{\link[permute]{how}} for details on permutation design
  specification, \code{\link[permute]{shuffleSet}} for the code used to
  generate a set of permutations, \code{\link[permute]{numPerms}} for
  a function to return the size of the set of possible permutations
  under the current design.
}

\references{

  Manly, B. F. J. (2006). \emph{Randomization, Bootstrap and Monte Carlo
  Methods in Biology}, Third Edition. Chapman and Hall/CRC.

  Phipson, B., & Smyth, G. K. (2010). Permutation P-values should never
  be zero: calculating exact P-values when permutations are randomly
  drawn. \emph{Statistical Applications in Genetics and Molecular
    Biology}, \strong{9}, Article 39. DOI: 10.2202/1544-6115.1585

  ter Braak, C. J. F. (1990). \emph{Update notes: CANOCO version
    3.1}. Wageningen: Agricultural Mathematics Group. (UR).

  See also:

  Davison, A. C., & Hinkley, D. V. (1997). \emph{Bootstrap Methods and
    their Application}. Cambridge University Press.
}
\author{ Gavin L. Simpson }
\keyword{multivariate}

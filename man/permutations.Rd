\name{permutations}
\alias{permutations}

\title{Permutation tests in Vegan}
\description{
  Unless stated otherwise, vegan currently provides for two types of
  permutation test:
  \enumerate{
    \item{Free permutation of \emph{DATA}, also known as randomisation,
      and}
    \item{Free permutation of \emph{DATA} within the levels of a factor
      variable.}
  }
  We use \emph{DATA} to mean either the observed data themselves or some
  function of the data, for example the residuals of an ordination model
  in the presence of covariables.
  
  The second type of permutation test above is available if the function
  providing the test accepts an argument \code{strata} or passes
  additional arguments (via \code{\dots}) to
  \code{\link{permuted.index}}.

  The Null hypothesis for these two types of permutation test assumes
  free exchangeability of \emph{DATA} (within the levels of
  \code{strata} if specified). Dependence between observations, such as
  that which arises due to spatial or temporal autocorrelation, or
  more-complicated experimental designs, such as split-plot designs,
  violates this fundamental assumption of the test and requires restricted
  permutation test designs. The next major version of Vegan will include
  infrastructure to handle these more complicated permutation designs.

  Again, unless otherwise stated in the help pages for specific
  functions, permutation tests in Vegan all follow the same
  format/structure:
  \enumerate{
    \item{An appropriate test statistic is chosen. Which statistic is
      chosen should be described on the help pages for individual
      functions.}
    \item{The value of the test statistic is evaluate for the observed
      data and analysis/model and recorded. Denote this value
      \eqn{x_0}{x[0]}.}
    \item{The \emph{DATA} are randomly permuted according to one of the
      above two schemes, and the value of the test statistic for this
      permutation is evaluated and recorded.}
    \item{Step 3 is repeated a total of \eqn{n} times, where \eqn{n} is
      the number of permutations requested. Denote these values as
      \eqn{x_i}{x[i]}, where \eqn{i = 1, ..., n}{{i = 1, \ldots, n}.}}
    \item{The values of the test statistic for the \eqn{n} permutations
      of the \emph{DATA} are added to the value of the test statistic
      for the observed data. These \emph{n + 1} values represent the
      \emph{Null} or \emph{randomisation} distribution of the test
      statistic. The observed value for the test statistic is included
      in the Null distribution because under the Null hypothesis being
      tested, the observed value is just a typical value of the test
      statistic, inherently no different from the values obtained via
      permutation of \emph{DATA}.}
    \item{The number of times that a value of the test statistic in the
      Null distribution is equal to or greater than the value of the
      test statistic for the observed data is recorded. Note the point
      mentioned in step 5 above; the Null distribution includes the
      \strong{observed} value of the test statistic. Denote this count
      as \eqn{N}.}
    \item{The permutation p-value is computed as
      \deqn{p = \frac{N}{n + 1}}{N / (n + 1)}}
  }
  The above description illustrates why the default number of
  permutations specified in Vegan functions takes values of 199 or 999
  for example. Once the observed value of the test statistic is added to
  this number of random permutations of \emph{DATA}, pretty p-values are
  achievable because \eqn{n + 1} becomes 200 or 1000, for example.

  The minimum achievable p-value is
  \deqn{p_{\mathrm{min}} = \frac{1}{n +1}}{p[min] = 1 / (n + 1)}
  
  A more common definition, in ecological circles, for \eqn{N} would be
  the number of \eqn{x_i}{x[i]} greater than or equal to
  \eqn{x_0}{x[0]}. The permutation p-value would then be defined as
  \deqn{p = \frac{N + 1}{n + 1}}{(N + 1) / (n + 1)}
  The + 1 in the numerator of the above equation represents the observed
  statistic \eqn{x_0}{x[0]}. The minimum p-value would then be defined as
  \deqn{p_{\mathrm{min}} = \frac{0 + 1}{n +1}}{p[min] = 0 + 1 / (n + 1)}
  However this definition discriminates between the observed
  statistic and the other \eqn{x_i}{x[i]}. Under the Null hypothesis
  there is no such distinction, hence we prefer the definintion used in
  the numbered steps above.

  One cannot simply increase the number of permutations
  (\eqn{n}) to achieve a potentially lower p-value unless the number of
  observations available permits such a number of permutations. This is
  unlikely to be a problem for all but the smallest data sets when
  free permutation (randomisation) is valid, but in designs where
  \code{strata} is specified and there are a low number of observations
  within each level of \code{strata}, there may not be as many actual
  permutations of the data as you might want.
  
  It is currently the responsibility of the user to determine the total
  number of possible permutations for their \emph{DATA}. No checks are
  made within Vegan functions to ensure a sensible number of
  permutations is chosen.

  Limits on the total number of permutations of \emph{DATA} are more
  severe in temporally or spatially ordered data or experimental designs
  with low replication. For example, a time series of \eqn{n = 100}
  observations has just 100 possible permutations \strong{including} the
  observed ordering.

  In situations where only a low number of permutations is possible due
  to the nature of \emph{DATA} or the experimental design, enumeration
  of all permutations becomes important and achievable
  computationally. Currently, Vegan does not include functions to
  perform complete enumeration of the set of possible
  permutations. The next major release of Vegan will include such
  functionality, however.
}

\seealso{
  \code{\link{permutest}}, \code{\link{permuted.index}}
}
%\references{
%}
\author{ Gavin Simpson }
\keyword{multivariate}

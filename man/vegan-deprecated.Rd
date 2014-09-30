\encoding{UTF-8}
\name{vegan-deprecated}
\alias{commsimulator}
\alias{density.adonis}
\alias{density.anosim}
\alias{density.mantel}
\alias{density.mrpp}
\alias{density.permutest.cca}
\alias{density.protest}
\alias{plot.vegandensity}
\alias{densityplot.adonis}
\alias{density.oecosimu}
\alias{densityplot.oecosimu}

\alias{vegan-deprecated}
%------ NOTE:  ../R/vegan-deprecated.R   must be synchronized with this!
\title{Deprecated Functions in vegan package}
%------ PLEASE: one \alias{.} for EACH ! (+ one \usage{} & \arguments{} for all)
\description{
  These functions are provided for compatibility with older versions of
  \pkg{vegan} only, and may be defunct as soon as the next release.
}
\usage{
commsimulator(x, method, thin=1)
\method{density}{adonis}(x, ...)
\method{plot}{vegandensity}(x, main = NULL, xlab = NULL, ylab = "Density", 
   type = "l", zero.line = TRUE, obs.line = TRUE, ...)
\method{densityplot}{adonis}(x, data, xlab = "Null", ...)
}

\arguments{
 \item{x}{Community data for \code{commsimulator},or an object to be
   handled by \code{density} or \code{densityplot}}
 ## commsimulator
 \item{method}{Null model method: either a name (character string) of
   a method defined in \code{\link{make.commsim}} or a
   \code{\link{commsim}} function.}
 \item{thin}{Number of discarded null communities between two
   evaluations of nestedness statistic in sequential methods
   \code{"swap"} and \code{"tswap"} (ignored with non-sequential
   methods)}
  ## density and densityplot
  \item{main, xlab, ylab, type, zero.line}{Arguments of
    \code{\link{plot.density}}, \code{\link[lattice]{densityplot}}.}
  \item{obs.line}{Draw vertical line for the observed
    statistic. Logical value \code{TRUE} draws a red line, and
    \code{FALSE} draws nothing. Alternatively, \code{obs.line} can be a
    definition of the colour used for the line, either as a numerical
    value from the \code{\link[grDevices]{palette}} or as the name of
    the colour, or other normal definition of the colour.}
  \item{data}{Ignored.}
  \item{\dots}{ Other arguments passed to functions. }

}

\details{
  %-- Either:
  %There are currently no deprecated functions in this package.
  %-- Or:
  %% explain *why* it's deprecated, and \code{\link{..}} to new
  
  Function \code{commsimulator} is replaced with
  \code{\link{make.commsim}} which defines the Null models, and
  functions \code{\link{nullmodel}} and
  \code{\link{simulate.nullmodel}} that check the input data and
  generate the Null model communities.  Function \code{commsimulator}
  was used to generate a single Null model for presence/absence
  (binary) data. Below is a copy of its original documentation in
  \code{\link{oecosimu}}, where it is now replaced with
  \code{\link{make.commsim}}, \code{\link{nullmodel}} and
  \code{\link{simulate.nullmodel}}.  Approximately the same
  documentation for these models is found in
  \code{\link{make.commsim}}.  (However, the random number sequences
  for model \code{r0} differ, and you must use \code{method = "r0_old"} 
  in \code{\link{make.commsim}} to reproduce the \code{commsimulator} 
  results.)

  Function \code{commsimulator} implements binary (presence/absence) 
  null models for community composition.
  The implemented models are \code{r00} which maintains the
  number of presences but fills these anywhere so that neither species
  (column) nor site (row) totals are preserved. Methods \code{r0},
  \code{r1} and \code{r2} maintain the site (row) frequencies. Method \code{r0}
  fills presences anywhere on the row with no respect to species (column)
  frequencies, \code{r1} uses column marginal 
  frequencies as probabilities, and \code{r2} uses squared column
  sums. Methods \code{r1} and \code{r2} try to simulate original species
  frequencies, but they are not strictly constrained. All these methods
  are reviewed by Wright et al. (1998). Method \code{c0} maintains
  species frequencies, but does not honour site (row) frequencies (Jonsson
  2001). 

  The other methods maintain both row and column frequencies.
  Methods \code{swap} and \code{tswap} implement sequential methods,
  where the matrix is changed only little in one step, but the changed
  matrix is used as an input if the next step.
  Methods \code{swap} and \code{tswap} inspect random 2x2 submatrices
  and if they are checkerboard units, the order of columns is
  swapped. This changes the matrix structure, but does not influence
  marginal sums (Gotelli & Entsminger
  2003). Method \code{swap} inspects submatrices so long that a swap
  can be done. \enc{Miklós}{Miklos} & Podani (2004) suggest that this may lead into
  biased sequences, since some columns or rows may be more easily
  swapped, and they suggest trying a fixed number of times and
  doing zero to many swaps at one step. This method is implemented by
  method \code{tswap} or trial swap. Function \code{commsimulator} makes
  only one trial swap in time (which probably does nothing),
  but \code{oecosimu} estimates how many
  submatrices are expected before finding a swappable checkerboard,
  and uses that ratio to thin the results, so that on average one swap
  will be found per step of \code{tswap}.  However, the checkerboard
  frequency probably changes during swaps, but this is not taken into
  account in estimating the \code{thin}.  One swap still changes the
  matrix only little, and it may be useful to 
  thin the results so that the statistic is only evaluated after
  \code{burnin} steps (and \code{thin}ned). 

  Methods \code{quasiswap} and \code{backtracking} are not sequential,
  but each call produces a matrix that is independent of previous
  matrices, and has the same marginal totals as the original data. The
  recommended method is \code{quasiswap} which is much faster because
  it is implemented in C. Method \code{backtracking} is provided for
  comparison, but it is so slow that it may be dropped from future
  releases of \pkg{vegan} (or also implemented in C).
  Method \code{quasiswap} (\enc{Miklós}{Miklos} & Podani 2004)
  implements a method where matrix is first filled 
  honouring row and column totals, but with integers that may be larger than
  one. Then the method inspects random 2x2 matrices and performs a
  quasiswap on them. Quasiswap is similar to ordinary swap, but it also
  can reduce numbers above one to ones maintaining marginal
  totals.
  Method \code{backtracking}
  implements a filling method with constraints both for row and column
  frequencies (Gotelli & Entsminger 2001). The matrix is first filled
  randomly using row and column frequencies as probabilities. Typically
  row and column sums are reached before all incidences are filled in.
  After that begins \dQuote{backtracking}, where some of the
  points are removed, and then filling is started again, and this
  backtracking is done so may times that all incidences will be filled
  into matrix. The \code{quasiswap} method is not sequential, but it produces
  a random incidence matrix with given marginal totals. 

  The \code{density} function can directly access permutation results
  of the same function as \code{permustats}.  The \code{density}
  function is identical to \code{\link{density.default}} and takes all
  its arguments, but adds the observed statistic to the result as item
  \code{"observed"}. The observed statistic is also put among the
  permuted values so that the results are consistent with significance
  tests. The \code{plot} method is similar to the default
  \code{\link{plot.density}}, but can also add the observed statistic
  to the graph as a vertical line.  In \code{\link{adonis}} it is also
  possible to use direclty \code{densityplot} function.

  The deprecated \code{density} and \code{densityplot} methods are
  replaced with similar methods for \code{\link{permustats}}. The
  \code{\link{permustats}} offers more powerful analysis tools for
  permutations, including \code{\link{summary.permustats}} giving
  \eqn{z} values (a.k.a. standardized effect sizes, SES), and Q-Q
  plots (\code{\link{qqnorm.permustats}},
  \code{\link{qqmath.permustats}}. Below the old documentation: 

  The density methods are available for \pkg{vegan} functions
  \code{\link{adonis}}, \code{\link{anosim}}, \code{\link{mantel}},
  \code{\link{mantel.partial}}, \code{\link{mrpp}},
  \code{\link{permutest.cca}}, and \code{\link{protest}}.  The
  \code{density} function for \code{\link{oecosimu}} is documented
  separately, and it is also used for \code{\link{adipart}},
  \code{\link{hiersimu}} and \code{\link{multipart}}.

  All \pkg{vegan} \code{density} functions return an object of class
  \code{"vegandensity"} inheriting from \code{\link{density}}, and can
  be plotted with its \code{plot} method.  This is identical to the
  standard \code{plot} of \code{densiy} objects, but can also add a
  vertical line for the observed statistic.

  Functions that can return several permuted statistics simultaneously
  also have \code{\link[lattice]{densityplot}} method
  (\code{\link{adonis}}, \code{\link{oecosimu}} and diversity
  partitioning functions based on \code{oecosimu}).  The standard
  \code{\link{density}} can only handle univariate data, and a warning
  is issued if the function is used for a model with several observed
  statistics

}

\references{
  Gotelli, N.J. & Entsminger, N.J. (2001). Swap and fill algorithms in
  null model analysis: rethinking the knight's tour. \emph{Oecologia}
  129, 281--291.

  Gotelli, N.J. & Entsminger, N.J. (2003). Swap algorithms in null model
  analysis. \emph{Ecology} 84, 532--535.

  Jonsson, B.G. (2001) A null model for randomization tests of
  nestedness in species assemblages. \emph{Oecologia} 127, 309--313.

  \enc{Miklós}{Miklos}, I. & Podani, J. (2004). Randomization of presence-absence
  matrices: comments and new algorithms. \emph{Ecology} 85, 86--92.

  Wright, D.H., Patterson, B.D., Mikkelson, G.M., Cutler, A. & Atmar,
  W. (1998). A comparative analysis of nested subset patterns of species
  composition. \emph{Oecologia} 113, 1--20.
  }

\seealso{
  \code{\link{Deprecated}}
}
\keyword{misc}

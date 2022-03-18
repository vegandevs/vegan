\encoding{UTF-8}
\name{commsim}
\alias{commsim}
\alias{make.commsim}
\alias{print.commsim}
\title{
Create an Object for Null Model Algorithms
}
\description{
The \code{commsim} function can be used to feed Null Model algorithms into
\code{\link{nullmodel}} analysis.
The \code{make.commsim} function returns various predefined algorithm types
(see Details).
These functions represent low level interface for community null model
infrastructure in \pkg{vegan} with the intent of extensibility,
and less emphasis on direct use by users.
}
\usage{
commsim(method, fun, binary, isSeq, mode)
make.commsim(method)
\method{print}{commsim}(x, ...)
}
\arguments{
  \item{method}{
Character, name of the algorithm.
}
  \item{fun}{
A function. For possible formal arguments of this function
see Details.
}
  \item{binary}{
Logical, if the algorithm applies to presence-absence or count matrices.
}
  \item{isSeq}{
Logical, if the algorithm is sequential (needs burnin and thinning) or not.
}
  \item{mode}{
Character, storage mode of the community matrix, either
\code{"integer"} or \code{"double"}.
}
  \item{x}{
An object of class \code{commsim}.
}
  \item{\dots}{
Additional arguments.
}
}
\details{
The function \code{fun} must return an array of \code{dim(nr, nc, n)},
and must take some of the following arguments:
\itemize{
  \item{\code{x}: }{input matrix,}
  \item{\code{n}: }{number of permuted matrices in output,}
  \item{\code{nr}: }{number of rows,}
  \item{\code{nc}: }{number of columns,}
  \item{\code{rs}: }{vector of row sums,}
  \item{\code{cs}: }{vector of column sums,}
  \item{\code{rf}: }{vector of row frequencies (non-zero cells),}
  \item{\code{cf}: }{vector of column frequencies (non-zero cells),}
  \item{\code{s}: }{total sum of \code{x},}
  \item{\code{fill}: }{matrix fill (non-zero cells),}
  \item{\code{thin}: }{thinning value for sequential algorithms,}
  \item{\code{...}: }{additional arguments.}
}

  You can define your own null model, but
  several null model algorithm are pre-defined and can be called by
  their name. The predefined algorithms are described in detail in the
  following chapters. The binary null models produce matrices of zeros
  (absences) and ones (presences) also when input matrix is
  quantitative. There are two types of quantitative data: Counts are
  integers with a natural unit so that individuals can be shuffled, but
  abundances can have real (floating point) values and do not have a
  natural subunit for shuffling. All quantitative models can handle
  counts, but only some are able to handle real values. Some of the null
  models are sequential so that the next matrix is derived from the
  current one. This makes models dependent from previous models, and usually
  you must thin these matrices and study the sequences for stability:
  see \code{oecosimu} for details and instructions.

  See Examples for structural constraints imposed by each algorithm and
  defining your own null model.

}
%% commsimulator

\section{Binary null models}{

  All binary null models preserve fill: number of presences or
  conversely the number of absences. The classic models may also
  preserve column (species) frequencies (\code{c0}) or row frequencies
  or species richness of each site (\code{r0}) and take into account
  commonness and rarity of species (\code{r1}, \code{r2}).  Algorithms
  \code{swap}, \code{tswap}, \code{curveball}, \code{quasiswap} and
  \code{backtracking} preserve both row and column frequencies. Three
  first ones are sequential but the two latter are non-sequential
  and produce independent matrices. Basic algorithms are reviewed by
  Wright et al. (1998).

\itemize{
  \item{\code{"r00"}: }{non-sequential algorithm for binary matrices
    that only  preserves the number of presences (fill).}

  \item{\code{"r0"}: }{non-sequential algorithm for binary
    matrices that preserves the site (row) frequencies.}

  \item{\code{"r1"}: }{non-sequential algorithm for binary matrices
    that preserves the site (row) frequencies, but uses column marginal
    frequencies as probabilities of selecting species.}

  \item{\code{"r2"}: }{non-sequential algorithm for binary matrices
    that preserves the site (row) frequencies, and uses squared column
    marginal frequencies as as probabilities of selecting species.}

  \item{\code{"c0"}: }{non-sequential algorithm for binary matrices
    that preserves species frequencies (Jonsson 2001). }

  \item{\code{"swap"}: }{sequential algorithm for binary matrices that
    changes the matrix structure, but does not influence marginal sums
    (Gotelli & Entsminger 2003).  This inspects \eqn{2 \times 2}{2 by
    2} submatrices so long that a swap can be done.}

  \item{\code{"tswap"}: }{sequential algorithm for binary matrices.
    Same as the \code{"swap"} algorithm, but it tries a fixed
    number of times and performs zero to many swaps at one step
    (according to the thin argument in the call). This
    approach was suggested by \enc{Miklós}{Miklos} & Podani (2004)
    because they found that ordinary swap may lead to biased
    sequences, since some columns or rows are more easily swapped.}

  \item{\code{"curveball"}: }{sequential method for binary matrices that
    implements the \sQuote{Curveball} algorithm of Strona et
    al. (2014). The algorithm selects two random rows and finds the set
    of unique species that occur only in one of these rows. The
    algorithm distributes the set of unique species to rows preserving
    the original row frequencies.  Zero to several species are swapped
    in one step, and usually the matrix is perturbed more strongly than
    in other sequential methods.}

  \item{\code{"quasiswap"}: }{non-sequential algorithm for binary
    matrices that implements a method where matrix is first filled
    honouring row and column totals, but with integers that may be
    larger than one.  Then the method inspects random
    \eqn{2 \times 2}{2 by 2} matrices and performs a quasiswap on
    them. In addition to ordinary swaps, quasiswap can reduce numbers
    above one to ones preserving marginal totals (\enc{Miklós}{Miklos} &
    Podani 2004). The method is non-sequential, but it accepts
    \code{thin} argument: the convergence is checked at every
    \code{thin} steps. This allows performing several ordinary swaps in
    addition to fill changing swaps which helps in reducing or removing
    the bias.}

  \item{\code{"greedyqswap"}: }{A greedy variant of quasiswap. In
    greedy step, one element of the \eqn{2 \times 2}{2 by 2} matrix is
    taken from \eqn{> 1} elements. The greedy steps are biased, but
    the method can be thinned, and only the first of \code{thin}
    steps is greedy. Even modest thinning (say \code{thin = 20})
    removes or reduces the bias, and \code{thin = 100} (1\% greedy
    steps) looks completely safe and still speeds up simulation. The
    code is experimental and it is provided here for further scrutiny,
    and should be tested for bias before use.}

  \item{\code{"backtracking"}: }{non-sequential algorithm for binary
    matrices that implements a filling method with constraints both for
    row and column frequencies (Gotelli & Entsminger 2001).  The matrix
    is first filled randomly, but typically row and column sums are
    reached before all incidences are filled in. After this begins
    "backtracking", where some of the incidences are removed, and
    filling is started again, and this backtracking is done so many
    times that all incidences will be filled into matrix.  The results
    may be biased and should be inspected carefully before use.}

}
}

\section{Quantitative Models for Counts with Fixed Marginal Sums}{

  These models shuffle individuals of counts and keep marginal sums
  fixed, but marginal frequencies are not preserved. Algorithm
  \code{r2dtable} uses standard \R function \code{\link{r2dtable}} also
  used for simulated \eqn{P}-values in \code{\link{chisq.test}}.
  Algorithm \code{quasiswap_count} uses the same, but preserves the
  original fill. Typically this means increasing numbers of zero cells
  and the result is zero-inflated with respect to \code{r2dtable}.

\itemize{

  \item{\code{"r2dtable"}: }{non-sequential algorithm for count
    matrices.  This algorithm keeps matrix sum and row/column sums
    constant. Based on \code{\link{r2dtable}}.}

  \item{\code{"quasiswap_count"}: }{non-sequential algorithm for count
    matrices.  This algorithm is similar as Carsten Dormann's
    \code{swap.web} function in the package
    \pkg{bipartite}. First, a random matrix is generated by the
    \code{\link{r2dtable}} function preserving row and column sums.  Then
    the original matrix fill is reconstructed by sequential steps to
    increase or decrease matrix fill in the random matrix. These steps
    are based on swapping \eqn{2 \times 2}{2 x 2} submatrices (see
    \code{"swap_count"} algorithm for details) to maintain row and
    column totals. }
}
}

\section{Quantitative Swap Models}{

  Quantitative swap models are similar to binary \code{swap}, but they
  swap the largest permissible value. The models in this section all
  maintain the fill and perform a quantitative swap only if this can
  be done without changing the fill. Single step of swap often changes
  the matrix very little. In particular, if cell counts are variable,
  high values change very slowly. Checking the chain stability and
  independence is even more crucial than in binary swap, and very
  strong \code{thin}ning is often needed. These models should never be
  used without inspecting their properties for the current data. These
  null models can also be defined using \code{\link{permatswap}}
  function.

 \itemize{

   \item{\code{"swap_count"}: }{sequential algorithm for count matrices.
    This algorithm find \eqn{2 \times 2}{2 x 2} submatrices that can be
    swapped leaving column and row totals and fill unchanged. The
    algorithm finds the largest value in the submatrix that can be
    swapped (\eqn{d}). Swap means that the values in diagonal or
    antidiagonal positions are decreased by \eqn{d}, while remaining
    cells are increased by \eqn{d}. A swap is made only if fill does not
    change.  }

   \item{\code{"abuswap_r"}: }{sequential algorithm for count or
    nonnegative real valued matrices with fixed row frequencies (see
    also \code{\link{permatswap}}).  The algorithm is similar to
    \code{swap_count}, but uses different swap value for each row of the
    \eqn{2 \times 2}{2 x 2} submatrix. Each step changes the the
    corresponding column sums, but honours matrix fill, row sums, and
    row/column frequencies (Hardy 2008; randomization scheme 2x).}

  \item{\code{"abuswap_c"}: }{sequential algorithm for count or
    nonnegative real valued matrices with fixed column frequencies
    (see also \code{\link{permatswap}}).  The algorithm is similar as
    the previous one, but operates on columns.  Each step changes the
    the corresponding row sums, but honours matrix fill, column sums,
    and row/column frequencies (Hardy 2008; randomization scheme 3x).}
    } }

\section{Quantitative Swap and Shuffle Models}{

  Quantitative Swap and Shuffle methods (\code{swsh} methods) preserve
  fill and column and row frequencies, and also either row or column
  sums. The methods first perform a binary \code{quasiswap} and then
  shuffle original quantitative data to non-zero cells. The
  \code{samp} methods shuffle original non-zero cell values and can be
  used also with non-integer data. The \code{both} methods
  redistribute individuals randomly among non-zero cells and can only
  be used with integer data. The shuffling is either free over the
  whole matrix, or within rows (\code{r} methods) or within columns
  (\code{c} methods). Shuffling within a row preserves row sums, and
  shuffling within a column preserves column sums. These models can
  also be defined with \code{\link{permatswap}}.

\itemize{

  \item{\code{"swsh_samp"}: }{non-sequential algorithm for
    quantitative data (either integer counts or non-integer values).
    Original non-zero values values are shuffled.}

  \item{\code{"swsh_both"}: }{non-sequential algorithm for count data.
    Individuals are shuffled freely over non-zero cells.}

  \item{\code{"swsh_samp_r"}: }{non-sequential algorithm for
    quantitative data.  Non-zero values (samples) are shuffled
    separately for each row.}

  \item{\code{"swsh_samp_c"}: }{non-sequential algorithm for
    quantitative data.  Non-zero values (samples) are shuffled
    separately for each column.}

  \item{\code{"swsh_both_r"}: }{non-sequential algorithm for count matrices.
    Individuals are shuffled freely for non-zero values within each row.}

  \item{\code{"swsh_both_c"}: }{non-sequential algorithm for count matrices.
    Individuals are shuffled freely for non-zero values with each column.}
}
}

\section{Quantitative Shuffle Methods}{

  Quantitative shuffle methods are generalizations of binary models
  \code{r00}, \code{r0} and \code{c0}.  The \code{_ind} methods
  shuffle individuals so that the grand sum, row sum or column sums
  are preserved.  These methods are similar as \code{r2dtable} but
  with still slacker constraints on marginal sums. The \code{_samp}
  and \code{_both} methods first apply the corresponding binary model
  with similar restriction on marginal frequencies and then distribute
  quantitative values over non-zero cells. The \code{_samp} models
  shuffle original cell values and can therefore handle also non-count
  real values. The \code{_both} models shuffle individuals among
  non-zero values. The shuffling is over the whole matrix in
  \code{r00_}, and within row in \code{r0_} and within column in
  \code{c0_} in all cases.

\itemize{
  \item{\code{"r00_ind"}: }{non-sequential algorithm for count matrices.
    This algorithm preserves grand sum and
    individuals are shuffled among cells of the matrix.}

  \item{\code{"r0_ind"}: }{non-sequential algorithm for count matrices.
    This algorithm preserves row sums and
    individuals are shuffled among cells of each row of the matrix.}

  \item{\code{"c0_ind"}: }{non-sequential algorithm for count matrices.
    This algorithm preserves column sums and
    individuals are shuffled among cells of each column of the matrix.}

  \item{\code{"r00_samp"}: }{non-sequential algorithm for count
    or nonnegative real valued (\code{mode = "double"}) matrices.
    This algorithm preserves grand sum and
    cells of the matrix are shuffled.}

  \item{\code{"r0_samp"}: }{non-sequential algorithm for count
    or nonnegative real valued (\code{mode = "double"}) matrices.
    This algorithm preserves row sums and
    cells within each row are shuffled.}

  \item{\code{"c0_samp"}: }{non-sequential algorithm for count
    or nonnegative real valued (\code{mode = "double"}) matrices.
    This algorithm preserves column sums constant and
    cells within each column are shuffled.}

  \item{\code{"r00_both"}: }{non-sequential algorithm for count matrices.
    This algorithm preserves grand sum and
    cells and individuals among cells of the matrix are shuffled.}

  \item{\code{"r0_both"}: }{non-sequential algorithm for count matrices.
    This algorithm preserves grand sum and
    cells and individuals among cells of each row are shuffled.}

  \item{\code{"c0_both"}: }{non-sequential algorithm for count matrices.
    This algorithm preserves grand sum and
    cells and individuals among cells of each column are shuffled.}
}
}

\value{
An object of class \code{commsim} with elements
corresponding to the arguments (\code{method}, \code{binary},
\code{isSeq}, \code{mode}, \code{fun}).

If the input of \code{make.comsimm} is a \code{commsim} object,
it is returned without further evaluation. If this is not the case,
the character \code{method} argument is matched against
predefined algorithm names. An error message is issued
if none such is found. If the \code{method} argument is missing,
the function returns names of all currently available
null model algorithms as a character vector.
}
\references{
  Gotelli, N.J. & Entsminger, N.J. (2001). Swap and fill algorithms in
  null model analysis: rethinking the knight's tour. \emph{Oecologia}
  129, 281--291.

  Gotelli, N.J. & Entsminger, N.J. (2003). Swap algorithms in null model
  analysis. \emph{Ecology} 84, 532--535.

  Hardy, O. J. (2008) Testing the spatial phylogenetic structure of
  local communities: statistical performances of different null models
  and test statistics on a locally neutral community.  \emph{Journal of
  Ecology} 96, 914--926.

  Jonsson, B.G. (2001) A null model for randomization tests of
  nestedness in species assemblages. \emph{Oecologia} 127, 309--313.

  \enc{Miklós}{Miklos}, I. & Podani, J. (2004). Randomization of
  presence-absence matrices: comments and new algorithms. \emph{Ecology}
  85, 86--92.

  Patefield, W. M. (1981) Algorithm AS159.  An efficient method of
  generating r x c tables with given row and column totals.
  \emph{Applied Statistics} 30, 91--97.

  Strona, G., Nappo, D., Boccacci, F., Fattorini, S. &
  San-Miguel-Ayanz, J. (2014). A fast and unbiased procedure to
  randomize ecological binary matrices with fixed row and column
  totals. \emph{Nature Communications} 5:4114
  \doi{10.1038/ncomms5114}.

  Wright, D.H., Patterson, B.D., Mikkelson, G.M., Cutler, A. & Atmar,
  W. (1998). A comparative analysis of nested subset patterns of species
  composition. \emph{Oecologia} 113, 1--20.
}

\author{
Jari Oksanen and Peter Solymos
}

\seealso{ See \code{\link{permatfull}}, \code{\link{permatswap}} for
alternative specification of quantitative null models. Function
\code{\link{oecosimu}} gives a higher-level interface for applying null
models in hypothesis testing and analysis of models. Function
\code{\link{nullmodel}} and \code{\link{simulate.nullmodel}} are used to
generate arrays of simulated null model matrices.  }

\examples{
## write the r00 algorithm
f <- function(x, n, ...)
    array(replicate(n, sample(x)), c(dim(x), n))
(cs <- commsim("r00", fun=f, binary=TRUE,
    isSeq=FALSE, mode="integer"))

## retrieving the sequential swap algorithm
(cs <- make.commsim("swap"))

## feeding a commsim object as argument
make.commsim(cs)

## making the missing c1 model using r1 as a template
##   non-sequential algorithm for binary matrices
##   that preserves the species (column) frequencies,
##   but uses row marginal frequencies
##   as probabilities of selecting sites
f <- function (x, n, nr, nc, rs, cs, ...) {
    out <- array(0L, c(nr, nc, n))
    J <- seq_len(nc)
    storage.mode(rs) <- "double"
    for (k in seq_len(n))
        for (j in J)
            out[sample.int(nr, cs[j], prob = rs), j, k] <- 1L
    out
}
cs <- make.commsim("r1")
cs$method <- "c1"
cs$fun <- f

## structural constraints
diagfun <- function(x, y) {
    c(sum = sum(y) == sum(x),
        fill = sum(y > 0) == sum(x > 0),
        rowSums = all(rowSums(y) == rowSums(x)),
        colSums = all(colSums(y) == colSums(x)),
        rowFreq = all(rowSums(y > 0) == rowSums(x > 0)),
        colFreq = all(colSums(y > 0) == colSums(x > 0)))
}
evalfun <- function(meth, x, n) {
    m <- nullmodel(x, meth)
    y <- simulate(m, nsim=n)
    out <- rowMeans(sapply(1:dim(y)[3],
        function(i) diagfun(attr(y, "data"), y[,,i])))
    z <- as.numeric(c(attr(y, "binary"), attr(y, "isSeq"),
        attr(y, "mode") == "double"))
    names(z) <- c("binary", "isSeq", "double")
    c(z, out)
}
x <- matrix(rbinom(10*12, 1, 0.5)*rpois(10*12, 3), 12, 10)
algos <- make.commsim()
a <- t(sapply(algos, evalfun, x=x, n=10))
print(as.table(ifelse(a==1,1,0)), zero.print = ".")
}
\keyword{ multivariate }
\keyword{ datagen }

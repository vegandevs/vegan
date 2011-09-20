\encoding{UTF-8}
\name{commsim}
\alias{commsim}
\alias{make.commsim}
\alias{print.commsim}
\title{
Create a Object for Null Model Algorithms
}
\description{
The \code{commsim} function can be used to feed Null Model algorithms into
\code{\link{nullmodel}} analysis.
The \code{make.commsim} function returns various predefined algorithm type
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
A function. For possible arguments see Details.
}
  \item{binary}{
Logical, if the algorithm applies to presence-abseence or count matrices.
}
  \item{isSeq}{
Logical, if the algorthm is sequential (needs burnin).
}
  \item{mode}{
Character, storage mode of community matrix, either 
\code{"integer"} or \code{"double"}.
}
  \item{x}{
A \code{commsim} object.
}
  \item{\dots}{
Additional arguments.
}
}
\details{
The function \code{fun} must take some of the
following arguments:
\itemize{
  item{\code{x}: }{input matrix,}
  item{\code{n}: }{number of permuted matrices in output,}
  item{\code{nr}: }{number of rows,}
  item{\code{nc}: }{number of columns,}
  item{\code{rs}: }{vector of row sums,}
  item{\code{cs}: }{vector of column sums,}
  item{\code{rf}: }{vector of row frequencies (nonzero cells),}
  item{\code{cf}: }{vector of column frequencies (nonzero cells),}
  item{\code{s}: }{total sum of \code{x}}
  item{\code{fill}: }{matrix fill,}
  item{\code{thin}: }{thinning value for sequential algorithms,}
  item{\code{...}: }{additional arguments.}
}
Autput must be an array of \code{dim(nr, nc, n)}.

The following algorithms are currently predefined:
\itemize{
%% commsimulator
  \item{\code{"r00"}: }{non-sequential algorithm for binary matrices described 
    in the page of \code{\link{commsimulator}}.
    This maintains the number of presences but fills 
    these anywhere so that neither species (column) nor 
    site (row) totals are preserved. See Wright et al. (1998) for review.}
  \item{\code{"r0"}: }{non-sequential algorithm for binary matrices described 
    in the page of \code{\link{commsimulator}}.
    This maintains the site (row) frequencies,
    fills presences anywhere on the row with no respect to 
    species (column) frequencies. See Wright et al. (1998) for review.}
  \item{\code{"r1"}: }{non-sequential algorithm for binary matrices described 
    in the page of \code{\link{commsimulator}}.
    This maintains the site (row) frequencies,
    uses column marginal frequencies as probabilities. 
    It tries to simulate original species frequencies, 
    but it is not strictly constrained. See Wright et al. (1998) for review.}
  \item{\code{"r2"}: }{non-sequential algorithm for binary matrices described 
    in the page of \code{\link{commsimulator}}.
    This maintains the site (row) frequencies,
    uses squared column sums as as probabilities. 
    It tries to simulate original species frequencies, 
    but it is not strictly constrained. See Wright et al. (1998) for review.}
  \item{\code{"c0"}: }{non-sequential algorithm for binary matrices described 
    in the page of \code{\link{commsimulator}}.
    This method maintains species frequencies, 
    but does not honour site (row) frequencies (Jonsson 2001). }
  \item{\code{"swap"}: }{sequential algorithm for binary matrices 
    described in the page of \code{\link{commsimulator}}.
    This changes the matrix structure, but does not influence 
    marginal sums (Gotelli & Entsminger 2003).
    This inspects 2 x 2 submatrices so long that a swap can be done.}
  \item{\code{"tswap"}: }{sequential algorithm for binary matrices 
    described in the page of \code{\link{commsimulator}}.
    Same as the \code{"swap"} algorithm, but it
    is trying a fixed number of times and doing zero to 
    many swaps at one step (according th thin argument in later call, 
    see Examples). This approach was suggested by
    \enc{Miklós}{Miklos} & Podani (2004) because they found that ordinary swap
    may lead into biased sequences, 
    since some columns or rows may be more easily swapped.}
  \item{\code{"quasiswap"}: }{non-sequential algorithm for binary matrices described 
    in the page of \code{\link{commsimulator}}.
    This method implements a 
    method where matrix is first filled honouring row and 
    column totals, but with integers that may be larger than one. 
    Then the method inspects random 2x2 matrices and performs a 
    quasiswap on them. Quasiswap is similar to ordinary swap, 
    but it also can reduce numbers above one to ones maintaining 
    marginal totals (\enc{Miklós}{Miklos} & Podani 2004). }
  \item{\code{"backtracking"}: }{non-sequential algorithm for binary matrices described 
    in the page of \code{\link{commsimulator}}.
    This method implements a filling method with constraints 
    both for row and column frequencies (Gotelli & Entsminger 2001). 
    The matrix is first filled randomly using row and column frequencies 
    as probabilities. Typically row and column sums are reached before 
    all incidences are filled in. After that begins "backtracking", 
    where some of the points are removed, and then filling is 
    started again, and this backtracking is done so may times that 
    all incidences will be filled into matrix.}
%% permatswap
  \item{\code{"swap_count"}: }{sequential algorithm for count matrices 
    described in the page of \code{\link{permatswap}}.
    This algorithm tries to find 2 x 2 submatrices 
    (identified by 2 random row and 2 random column indices), 
    that can be swapped in order to leave column and row totals 
    and fill unchanged. First, the algorithm finds the largest 
    value in the submatrix that can be swapped (\eqn{d}) 
    and whether in diagonal or antidiagonal way. 
    Submatrices that contain values larger than zero in either 
    diagonal or antidiagonal position can be swapped. 
    Swap means that the values in diagonal or antidiagonal 
    positions are decreased by \eqn{d}, while remaining cells 
    are increased by \eqn{d}. A swap is made only if fill doesn't change. 
    \bold{WARNING}: according to simulations, 
    this algorithm seems to be biased and non random, 
    thus its use should be avoided!}
  \item{\code{"quasiswap_count"}: }{non-sequential algorithm for count matrices 
    described in the page of \code{\link{permatswap}}.
    This algorithm uses the same trick as 
    Carsten Dormann's \code{\link[bipartite]{swap.web}} 
    function in the package \pkg{bipartite}. First, a 
    random matrix is generated by the \code{\link{r2dtable}} 
    function retaining row and column sums. 
    Then the original matrix fill is reconstructed by 
    sequential steps to increase or decrease matrix fill in 
    the random matrix. These steps are based on swapping 2 x 2 
    submatrices (see \code{"swap_count"} algorithm for details) 
    to maintain row and column totals. }
  \item{\code{"swsh_samp"}: }{non-sequential algorithm for count matrices 
    described in the page of \code{\link{permatswap}}.
    The algorithm is a hybrid algorithm. 
    First, it makes binary quasiswaps to keep row and column 
    incidences constant, then non-zero values are 
    modified by shuffling original nonzero samples among the new
    nonzero cells of the binary matrix.}
  \item{\code{"swsh_both"}: }{non-sequential algorithm for count matrices 
    described in the page of \code{\link{permatswap}}.
    The algorithm is a hybrid algorithm. 
    First, it makes binary quasiswaps to keep row and column 
    incidences constant, then non-zero values are 
    modified by shuffling original nonzero samples among the new
    nonzero cells of the binary matrix, individuals of the
    nonzero cells are also shuffled.}
  \item{\code{"swsh_samp_r"}: }{non-sequential algorithm for count matrices 
    described in the page of \code{\link{permatswap}}.
    The algorithm is a hybrid algorithm. 
    First, it makes binary quasiswaps to keep row and column 
    incidences constant, then non-zero values are 
    modified by shuffling original nonzero samples among the new
    nonzero cells of the binary matrix, separately for each rows.}
  \item{\code{"swsh_samp_c"}: }{non-sequential algorithm for count matrices 
    described in the page of \code{\link{permatswap}}.
    The algorithm is a hybrid algorithm. 
    First, it makes binary quasiswaps to keep row and column 
    incidences constant, then non-zero values are 
    modified by shuffling original nonzero samples among the new
    nonzero cells of the binary matrix, separately for each columns.}
  \item{\code{"swsh_both_r"}: }{non-sequential algorithm for count matrices 
    described in the page of \code{\link{permatswap}}.
    The algorithm is a hybrid algorithm. 
    First, it makes binary quasiswaps to keep row and column 
    incidences constant, then non-zero values are 
    modified by shuffling original nonzero samples among the new
    nonzero cells of the binary matrix, individuals of the
    nonzero cells are also shuffled, separately for each rows.}
  \item{\code{"swsh_both_c"}: }{non-sequential algorithm for count matrices 
    described in the page of \code{\link{permatswap}}.
    The algorithm is a hybrid algorithm. 
    First, it makes binary quasiswaps to keep row and column 
    incidences constant, then non-zero values are 
    modified by shuffling original nonzero samples among the new
    nonzero cells of the binary matrix, individuals of the
    nonzero cells are also shuffled, separately for each columns.}
  \item{\code{"abuswap_r"}: }{sequential algorithm for count matrices 
    described in the page of \code{\link{permatswap}} (\code{"abuswap"} 
    algorithm with fixed row freqencies).
    The algorithm produces matrices by swapping
    2 x 2 submatrices. Each step changes the the corresponding
    column sums, but honours matrix fill, row sums,
    and row/column frequencies, 
    as described in Hardy (2008; randomization scheme 2x).}
  \item{\code{"abuswap_c"}: }{sequential algorithm for count matrices 
    described in the page of \code{\link{permatswap}} (\code{"abuswap"} 
    algorithm with fixed column freqencies).
    The algorithm produces matrices by swapping
    2 x 2 submatrices. Each step changes the the corresponding
    row sums, but honours matrix fill, column sums,
    and row/column frequencies, 
    as described in Hardy (2008; randomization scheme 3x).}
%% permatfull
  \item{\code{"r2dtable"}: }{non-sequential algorithm for count matrices 
    described in the page of \code{\link{r2dtable}} and \code{\link{permatfull}}.
    This algorithm keeps matrix sum and row/column sums constant.}
  \item{\code{"r00_ind"}: }{non-sequential algorithm for count matrices 
    described in the page of \code{\link{permatfull}}.
    This algorithm keeps total sum constant,
    individuals are shullfed among cells of the matrix.}
  \item{\code{"r0_ind"}: }{non-sequential algorithm for count matrices 
    described in the page of \code{\link{permatfull}}.
    This algorithm keeps row sums constant,
    individuals are shullfed among cells of each row of the matrix.}
  \item{\code{"c0_ind"}: }{non-sequential algorithm for count matrices 
    described in the page of \code{\link{permatfull}}.
    This algorithm keeps column sums constant,
    individuals are shullfed among cells of each column of the matrix.}
  \item{\code{"r00_samp"}: }{non-sequential algorithm for count matrices 
    described in the page of \code{\link{permatfull}}.
    This algorithm keeps total sum constant,
    cells of the matrix are shuffled.}
  \item{\code{"r0_samp"}: }{non-sequential algorithm for count matrices 
    described in the page of \code{\link{permatfull}}.
    This algorithm keeps row sums constant,
    cells within each row are shuffled.}
  \item{\code{"c0_samp"}: }{non-sequential algorithm for count matrices 
    described in the page of \code{\link{permatfull}}.
    This algorithm keeps column sums constant,
    cells within each columns are shuffled.}
  \item{\code{"r00_both"}: }{non-sequential algorithm for count matrices 
    described in the page of \code{\link{permatfull}}.
    This algorithm keeps total sum constant,
    cells and individuals among cells of the matrix are shuffled.}
  \item{\code{"r0_both"}: }{non-sequential algorithm for count matrices 
    described in the page of \code{\link{permatfull}}.
    This algorithm keeps total sum constant,
    cells and individuals among cells of each row are shuffled.}
  \item{\code{"c0_both"}: }{non-sequential algorithm for count matrices 
    described in the page of \code{\link{permatfull}}.
    This algorithm keeps total sum constant,
    cells and individuals among cells of each column are shuffled.}
}
For structural constraints imposed by each algorithm, see Examples.
}
\value{
An object of class \code{"commsim"} with elements 
corresponding to the arguments (\code{method}, \code{binary}, 
\code{isSeq}, \code{fun}).

IF the input of \code{make.comsimm} is a \code{commsim} object,
it is returned without further evaluation. If this is not the case,
the character \code{method} argument is matched agains
predefined algorithm names. An error message is issued
if none such is found.
}
\references{
  Hardy, O. J. (2008) 
  Testing the spatial phylogenetic structure of local communities: 
  statistical performances of different null models 
  and test statistics on a locally neutral community. 
  \emph{Journal of Ecology} 96, 914--926.

  Patefield, W. M. (1981) Algorithm AS159. 
  An efficient method of generating r x c tables with 
  given row and column totals. 
  \emph{Applied Statistics} 30, 91--97.

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
\author{
Peter Solymos \email{solymos@ualberta.ca}
}
\seealso{
\code{\link{commsimulator}}, \code{\link{permatfull}}, 
\code{\link{permatswap}}, \code{\link{oecosimu}}.
}
\examples{
## write the r00 algorithm
f <- function(x, n, ...) 
    array(replicate(n, sample(x)), c(dim(x), n))
(cs <- commsim("r00", fun=f, binary=TRUE, 
    isSeq=FALSE, mode="integer"))

## retrieving the sequential swap algorithm
make.commsim("swap")

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
    z <- as.numeric(c(attr(y, "binary"), attr(y, "isSeq")))
    names(z) <- c("binary", "isSeq")
    c(z, out)
}
x <- matrix(rbinom(10*12, 1, 0.5)*rpois(10*12, 3), 12, 10)
a <- t(sapply(c("r00","r0","r1","r2","c0",
    "swap","tswap","quasiswap","backtrack",
    "r2dtable","swap_count","quasiswap_count",
    "swsh_samp","swsh_both","abuswap_r","abuswap_c",
    "swsh_samp_r","swsh_samp_c","swsh_both_r","swsh_both_c",
    "r00_ind","r0_ind","c0_ind",
    "r00_samp","r0_samp","c0_samp",
    "r00_both","r0_both","c0_both"), 
    evalfun, x=x, n=10))
print(as.table(ifelse(a==1,1,0)), zero.print = ".")
}
\keyword{ multivariate }
\keyword{ datagen }

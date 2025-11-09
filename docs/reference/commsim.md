# Create an Object for Null Model Algorithms

The `commsim` function can be used to feed Null Model algorithms into
[`nullmodel`](https://vegandevs.github.io/vegan/reference/nullmodel.md)
analysis. The `make.commsim` function returns various predefined
algorithm types (see Details). These functions represent low level
interface for community null model infrastructure in vegan with the
intent of extensibility, and less emphasis on direct use by users.

## Usage

``` r
commsim(method, fun, binary, isSeq, mode)
make.commsim(method)
# S3 method for class 'commsim'
print(x, ...)
```

## Arguments

- method:

  Character, name of the algorithm.

- fun:

  A function. For possible formal arguments of this function see
  Details.

- binary:

  Logical, if the algorithm applies to presence-absence or count
  matrices.

- isSeq:

  Logical, if the algorithm is sequential (needs burnin and thinning) or
  not.

- mode:

  Character, storage mode of the community matrix, either `"integer"` or
  `"double"`.

- x:

  An object of class `commsim`.

- ...:

  Additional arguments.

## Details

The function `fun` must return an array of `dim(nr, nc, n)`, and must
take some of the following arguments:

- `x`: :

  input matrix,

- `n`: :

  number of permuted matrices in output,

- `nr`: :

  number of rows,

- `nc`: :

  number of columns,

- `rs`: :

  vector of row sums,

- `cs`: :

  vector of column sums,

- `rf`: :

  vector of row frequencies (non-zero cells),

- `cf`: :

  vector of column frequencies (non-zero cells),

- `s`: :

  total sum of `x`,

- `fill`: :

  matrix fill (non-zero cells),

- `thin`: :

  thinning value for sequential algorithms,

- `...`: :

  additional arguments.

You can define your own null model, but several null model algorithm are
pre-defined and can be called by their name. The predefined algorithms
are described in detail in the following chapters. The binary null
models produce matrices of zeros (absences) and ones (presences) also
when input matrix is quantitative. There are two types of quantitative
data: Counts are integers with a natural unit so that individuals can be
shuffled, but abundances can have real (floating point) values and do
not have a natural subunit for shuffling. All quantitative models can
handle counts (integers), but only some are able to handle real values.
Some of the null models are sequential so that the next matrix is
derived from the current one. This makes models dependent from previous
models, and usually you must thin these matrices and study the sequences
for stability: see
[`oecosimu`](https://vegandevs.github.io/vegan/reference/oecosimu.md)
for details and instructions.

See Examples for structural constraints imposed by each algorithm and
defining your own null model.

## Binary null models

All binary null models preserve fill: number of presences or conversely
the number of absences. The classic models may also preserve column
(species) frequencies (`c0`) or row frequencies or species richness of
each site (`r0`) and take into account commonness and rarity of species
(`r1`, `r2`). Algorithms `swap`, `tswap`, `curveball`, `quasiswap` and
`backtrack` preserve both row and column frequencies. Three first ones
are sequential but the two latter are non-sequential and produce
independent matrices. Basic algorithms are reviewed by Wright et al.
(1998).

- `"r00"`: :

  non-sequential algorithm for binary matrices that only preserves the
  number of presences (fill).

- `"r0"`: :

  non-sequential algorithm for binary matrices that preserves the site
  (row) frequencies.

- `"r1"`: :

  non-sequential algorithm for binary matrices that preserves the site
  (row) frequencies, but uses column marginal frequencies as
  probabilities of selecting species.

- `"r2"`: :

  non-sequential algorithm for binary matrices that preserves the site
  (row) frequencies, and uses squared column marginal frequencies as as
  probabilities of selecting species.

- `"c0"`: :

  non-sequential algorithm for binary matrices that preserves species
  frequencies (Jonsson 2001).

- `"swap"`: :

  sequential algorithm for binary matrices that changes the matrix
  structure, but does not influence marginal sums (Gotelli & Entsminger
  2003). This inspects \\2 \times 2\\ submatrices so long that a swap
  can be done. Some arbitrary matrix with (nearly) complete fill cannot
  be swapped; function
  [`nestedchecker`](https://vegandevs.github.io/vegan/reference/nestedtemp.md)
  gives the number of swappable submatrices. This also applies to other
  `swap` models.

- `"tswap"`: :

  sequential algorithm for binary matrices. Same as the `"swap"`
  algorithm, but it tries a fixed number of times and performs zero to
  many swaps at one step (according to the thin argument in the call).
  This approach was suggested by Miklós & Podani (2004) because they
  found that ordinary swap may lead to biased sequences, since some
  columns or rows are more easily swapped.

- `"curveball"`: :

  sequential method for binary matrices that implements the ‘Curveball’
  algorithm of Strona et al. (2014). The algorithm selects two random
  rows and finds the set of unique species that occur only in one of
  these rows. The algorithm distributes the set of unique species to
  rows preserving the original row frequencies. Zero to several species
  are swapped in one step, and usually the matrix is perturbed more
  strongly than in other sequential methods.

- `"quasiswap"`: :

  non-sequential algorithm for binary matrices that implements a method
  where matrix is first filled honouring row and column totals, but with
  integers that may be larger than one. Then the method inspects random
  \\2 \times 2\\ matrices and performs a quasiswap on them. In addition
  to ordinary swaps, quasiswap can reduce numbers above one to ones
  preserving marginal totals (Miklós & Podani 2004). The method is
  non-sequential, but it accepts `thin` argument: the convergence is
  checked at every `thin` steps. This allows performing several ordinary
  swaps in addition to fill changing swaps which helps in reducing or
  removing the bias.

- `"greedyqswap"`: :

  A greedy variant of quasiswap. In greedy step, one element of the \\2
  \times 2\\ matrix is taken from \\\> 1\\ elements. The greedy steps
  are biased, but the method can be thinned, and only the first of
  `thin` steps is greedy. Even modest thinning (say `thin = 20`) removes
  or reduces the bias, and `thin = 100` (1% greedy steps) looks
  completely safe and still speeds up simulation. The code is
  experimental and it is provided here for further scrutiny, and should
  be tested for bias before use.

- `"backtrack"`: :

  non-sequential algorithm for binary matrices that implements a filling
  method with constraints both for row and column frequencies (Gotelli &
  Entsminger 2001). The matrix is first filled randomly, but typically
  row and column sums are reached before all incidences are filled in.
  After this the function “backtracks” removing some incidences and
  starting to fill again. This backtracking is done so many times that
  all incidences will be filled into matrix. The results may be biased
  and should be inspected carefully before use.

## Quantitative Models for Counts with Fixed Marginal Sums

These models shuffle individuals of counts and keep marginal sums fixed,
but marginal frequencies are not preserved. Algorithm `r2dtable` uses
standard R function [`r2dtable`](https://rdrr.io/r/stats/r2dtable.html)
also used for simulated \\P\\-values in
[`chisq.test`](https://rdrr.io/r/stats/chisq.test.html). Algorithm
`quasiswap_count` uses the same, but preserves the original fill.
Typically this means increasing numbers of zero cells and the result is
zero-inflated with respect to `r2dtable`.

- `"r2dtable"`: :

  non-sequential algorithm for count matrices. This algorithm keeps
  matrix sum and row/column sums constant. Based on
  [`r2dtable`](https://rdrr.io/r/stats/r2dtable.html).

- `"quasiswap_count"`: :

  non-sequential algorithm for count matrices. This algorithm is similar
  as Carsten Dormann's `swap.web` function in the package bipartite.
  First, a random matrix is generated by the
  [`r2dtable`](https://rdrr.io/r/stats/r2dtable.html) function
  preserving row and column sums. Then the original matrix fill is
  reconstructed by sequential steps to increase or decrease matrix fill
  in the random matrix. These steps are based on swapping \\2 \times 2\\
  submatrices (see `"swap_count"` algorithm for details) to maintain row
  and column totals.

## Quantitative Swap Models

Quantitative swap models are similar to binary `swap`, but they swap the
largest permissible value. The models in this section all maintain the
fill and perform a quantitative swap only if this can be done without
changing the fill. Single step of swap often changes the matrix very
little. In particular, if cell counts are variable, high values change
very slowly. Checking the chain stability and independence is even more
crucial than in binary swap, and very strong `thin`ning is often needed.
These models should never be used without inspecting their properties
for the current data. These null models can also be defined using
[`permatswap`](https://vegandevs.github.io/vegan/reference/permatfull.md)
function.

- `"swap_count"`: :

  sequential algorithm for count matrices. This algorithm find \\2
  \times 2\\ submatrices that can be swapped leaving column and row
  totals and fill unchanged. The algorithm finds the largest value in
  the submatrix that can be swapped (\\d\\). Swap means that the values
  in diagonal or antidiagonal positions are decreased by \\d\\, while
  remaining cells are increased by \\d\\. A swap is made only if fill
  does not change.

- `"abuswap_r"`: :

  sequential algorithm for count or nonnegative real valued matrices
  with fixed row frequencies (see also
  [`permatswap`](https://vegandevs.github.io/vegan/reference/permatfull.md)).
  The algorithm is similar to `swap_count`, but uses different swap
  value for each row of the \\2 \times 2\\ submatrix. Each step changes
  the the corresponding column sums, but honours matrix fill, row sums,
  and row/column frequencies (Hardy 2008; randomization scheme 2x).

- `"abuswap_c"`: :

  sequential algorithm for count or nonnegative real valued matrices
  with fixed column frequencies (see also
  [`permatswap`](https://vegandevs.github.io/vegan/reference/permatfull.md)).
  The algorithm is similar as the previous one, but operates on columns.
  Each step changes the the corresponding row sums, but honours matrix
  fill, column sums, and row/column frequencies (Hardy 2008;
  randomization scheme 3x).

## Quantitative Swap and Shuffle Models

Quantitative Swap and Shuffle methods (`swsh` methods) preserve fill and
column and row frequencies, and also either row or column sums. The
methods first perform a binary `quasiswap` and then shuffle original
quantitative data to non-zero cells. The `samp` methods shuffle original
non-zero cell values and can be used also with non-integer data. The
`both` methods redistribute individuals randomly among non-zero cells
and can only be used with integer data. The shuffling is either free
over the whole matrix, or within rows (`r` methods) or within columns
(`c` methods). Shuffling within a row preserves row sums, and shuffling
within a column preserves column sums. These models can also be defined
with
[`permatswap`](https://vegandevs.github.io/vegan/reference/permatfull.md).

- `"swsh_samp"`: :

  non-sequential algorithm for quantitative data (either integer counts
  or non-integer values). Original non-zero values values are shuffled.

- `"swsh_both"`: :

  non-sequential algorithm for count data. Individuals are shuffled
  freely over non-zero cells.

- `"swsh_samp_r"`: :

  non-sequential algorithm for quantitative data. Non-zero values
  (samples) are shuffled separately for each row.

- `"swsh_samp_c"`: :

  non-sequential algorithm for quantitative data. Non-zero values
  (samples) are shuffled separately for each column.

- `"swsh_both_r"`: :

  non-sequential algorithm for count matrices. Individuals are shuffled
  freely for non-zero values within each row.

- `"swsh_both_c"`: :

  non-sequential algorithm for count matrices. Individuals are shuffled
  freely for non-zero values with each column.

## Quantitative Shuffle Methods

Quantitative shuffle methods are generalizations of binary models `r00`,
`r0` and `c0`. The `_ind` methods shuffle individuals so that the grand
sum, row sum or column sums are preserved. These methods are similar as
`r2dtable` but with still slacker constraints on marginal sums. The
`_samp` and `_both` methods first apply the corresponding binary model
with similar restriction on marginal frequencies and then distribute
quantitative values over non-zero cells. The `_samp` models shuffle
original cell values and can therefore handle also non-count real
values. The `_both` models shuffle individuals among non-zero values.
The shuffling is over the whole matrix in `r00_`, and within row in
`r0_` and within column in `c0_` in all cases.

- `"r00_ind"`: :

  non-sequential algorithm for count matrices. This algorithm preserves
  grand sum and individuals are shuffled among cells of the matrix.

- `"r0_ind"`: :

  non-sequential algorithm for count matrices. This algorithm preserves
  row sums and individuals are shuffled among cells of each row of the
  matrix.

- `"c0_ind"`: :

  non-sequential algorithm for count matrices. This algorithm preserves
  column sums and individuals are shuffled among cells of each column of
  the matrix.

- `"r00_samp"`: :

  non-sequential algorithm for count or nonnegative real valued
  (`mode = "double"`) matrices. This algorithm preserves grand sum and
  cells of the matrix are shuffled.

- `"r0_samp"`: :

  non-sequential algorithm for count or nonnegative real valued
  (`mode = "double"`) matrices. This algorithm preserves row sums and
  cells within each row are shuffled.

- `"c0_samp"`: :

  non-sequential algorithm for count or nonnegative real valued
  (`mode = "double"`) matrices. This algorithm preserves column sums
  constant and cells within each column are shuffled.

- `"r00_both"`: :

  non-sequential algorithm for count matrices. This algorithm preserves
  grand sum and cells and individuals among cells of the matrix are
  shuffled.

- `"r0_both"`: :

  non-sequential algorithm for count matrices. This algorithm preserves
  grand sum and cells and individuals among cells of each row are
  shuffled.

- `"c0_both"`: :

  non-sequential algorithm for count matrices. This algorithm preserves
  grand sum and cells and individuals among cells of each column are
  shuffled.

## Value

An object of class `commsim` with elements corresponding to the
arguments (`method`, `binary`, `isSeq`, `mode`, `fun`).

If the input of `make.comsimm` is a `commsim` object, it is returned
without further evaluation. If this is not the case, the character
`method` argument is matched against predefined algorithm names. An
error message is issued if none such is found. If the `method` argument
is missing, the function returns names of all currently available null
model algorithms as a character vector.

## References

Gotelli, N.J. & Entsminger, N.J. (2001). Swap and fill algorithms in
null model analysis: rethinking the knight's tour. *Oecologia* 129,
281–291.

Gotelli, N.J. & Entsminger, N.J. (2003). Swap algorithms in null model
analysis. *Ecology* 84, 532–535.

Hardy, O. J. (2008) Testing the spatial phylogenetic structure of local
communities: statistical performances of different null models and test
statistics on a locally neutral community. *Journal of Ecology* 96,
914–926.

Jonsson, B.G. (2001) A null model for randomization tests of nestedness
in species assemblages. *Oecologia* 127, 309–313.

Miklós, I. & Podani, J. (2004). Randomization of presence-absence
matrices: comments and new algorithms. *Ecology* 85, 86–92.

Patefield, W. M. (1981) Algorithm AS159. An efficient method of
generating r x c tables with given row and column totals. *Applied
Statistics* 30, 91–97.

Strona, G., Nappo, D., Boccacci, F., Fattorini, S. & San-Miguel-Ayanz,
J. (2014). A fast and unbiased procedure to randomize ecological binary
matrices with fixed row and column totals. *Nature Communications*
5:4114 [doi:10.1038/ncomms5114](https://doi.org/10.1038/ncomms5114) .

Wright, D.H., Patterson, B.D., Mikkelson, G.M., Cutler, A. & Atmar, W.
(1998). A comparative analysis of nested subset patterns of species
composition. *Oecologia* 113, 1–20.

## Author

Jari Oksanen and Peter Solymos

## See also

See
[`permatfull`](https://vegandevs.github.io/vegan/reference/permatfull.md),
[`permatswap`](https://vegandevs.github.io/vegan/reference/permatfull.md)
for alternative specification of quantitative null models. Function
[`oecosimu`](https://vegandevs.github.io/vegan/reference/oecosimu.md)
gives a higher-level interface for applying null models in hypothesis
testing and analysis of models. Function
[`nullmodel`](https://vegandevs.github.io/vegan/reference/nullmodel.md)
and
[`simulate.nullmodel`](https://vegandevs.github.io/vegan/reference/nullmodel.md)
are used to generate arrays of simulated null model matrices.

## Examples

``` r
## write the r00 algorithm
f <- function(x, n, ...)
    array(replicate(n, sample(x)), c(dim(x), n))
(cs <- commsim("r00", fun=f, binary=TRUE,
    isSeq=FALSE, mode="integer"))
#> An object of class “commsim” 
#> ‘r00’ method (binary, non-sequential, integer mode)
#> 

## retrieving the sequential swap algorithm
(cs <- make.commsim("swap"))
#> An object of class “commsim” 
#> ‘swap’ method (binary, sequential, integer mode)
#> 

## feeding a commsim object as argument
make.commsim(cs)
#> An object of class “commsim” 
#> ‘swap’ method (binary, sequential, integer mode)
#> 

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
#>                 binary isSeq double sum fill rowSums colSums rowFreq colFreq
#> r00                  1     .      .   1    1       .       .       .       .
#> c0                   1     .      .   1    1       .       1       .       1
#> r0                   1     .      .   1    1       1       .       1       .
#> r1                   1     .      .   1    1       1       .       1       .
#> r2                   1     .      .   1    1       1       .       1       .
#> quasiswap            1     .      .   1    1       1       1       1       1
#> greedyqswap          1     .      .   1    1       1       1       1       1
#> swap                 1     1      .   1    1       1       1       1       1
#> tswap                1     1      .   1    1       1       1       1       1
#> curveball            1     1      .   1    1       1       1       1       1
#> backtrack            1     .      .   1    1       1       1       1       1
#> r2dtable             .     .      .   1    .       1       1       .       .
#> swap_count           .     1      .   1    1       1       1       .       .
#> quasiswap_count      .     .      .   1    1       1       1       .       .
#> swsh_samp            .     .      1   1    1       .       .       1       1
#> swsh_both            .     .      .   1    1       .       .       1       1
#> swsh_samp_r          .     .      1   1    1       1       .       1       1
#> swsh_samp_c          .     .      1   1    1       .       1       1       1
#> swsh_both_r          .     .      .   1    1       1       .       1       1
#> swsh_both_c          .     .      .   1    1       .       1       1       1
#> abuswap_r            .     1      1   1    1       1       .       1       1
#> abuswap_c            .     1      1   1    1       .       1       1       1
#> r00_samp             .     .      1   1    1       .       .       .       .
#> c0_samp              .     .      1   1    1       .       1       .       1
#> r0_samp              .     .      1   1    1       1       .       1       .
#> r00_ind              .     .      .   1    .       .       .       .       .
#> c0_ind               .     .      .   1    .       .       1       .       .
#> r0_ind               .     .      .   1    .       1       .       .       .
#> r00_both             .     .      .   1    1       .       .       .       .
#> c0_both              .     .      .   1    1       .       1       .       1
#> r0_both              .     .      .   1    1       1       .       1       .
```

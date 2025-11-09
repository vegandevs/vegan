# Null Model and Simulation

The `nullmodel` function creates an object which can serve as a basis
for Null Model simulation via the
[`simulate`](https://rdrr.io/r/stats/simulate.html) method. The
[`update`](https://rdrr.io/r/stats/update.html) method updates the
nullmodel object without sampling (effective for sequential algorithms).
`smbind` binds together multiple `simmat` objects.

## Usage

``` r
nullmodel(x, method)
# S3 method for class 'nullmodel'
print(x, ...)
# S3 method for class 'nullmodel'
simulate(object, nsim = 1, seed = NULL,
    burnin = 0, thin = 1, ...)
# S3 method for class 'nullmodel'
update(object, nsim = 1, seed = NULL, ...)
# S3 method for class 'simmat'
print(x, ...)
smbind(object, ..., MARGIN, strict = TRUE)
```

## Arguments

- x:

  A community matrix. For the `print` method, it is an object to be
  printed.

- method:

  Character, specifying one of the null model algorithms listed on the
  help page of
  [`commsim`](https://vegandevs.github.io/vegan/reference/commsim.md).
  It can be a user supplied object of class `commsim`.

- object:

  An object of class `nullmodel` returned by the function `nullmodel`.
  In case of `smbind` it is a `simmat` object as returned by the
  `update` or `simulate` methods.

- nsim:

  Positive integer, the number of simulated matrices to return. For the
  `update` method, it is the number of burnin steps made for sequential
  algorithms to update the status of the input model `object`.

- seed:

  An object specifying if and how the random number generator should be
  initialized ("seeded"). Either `NULL` or an integer that will be used
  in a call to [`set.seed`](https://rdrr.io/r/base/Random.html) before
  simulating the matrices. If set, the value is saved as the `"seed"`
  attribute of the returned value. The default, `NULL` will not change
  the random generator state, and return
  [`.Random.seed`](https://rdrr.io/r/base/Random.html) as the `"seed"`
  attribute, see Value.

- burnin:

  Nonnegative integer, specifying the number of steps discarded before
  starting simulation. Active only for sequential null model algorithms.
  Ignored for non-sequential null model algorithms.

- thin:

  Positive integer, number of simulation steps made between each
  returned matrix. Active only for sequential null model algorithms.
  Ignored for non-sequential null model algorithms.

- MARGIN:

  Integer, indicating the dimension over which multiple `simmat` objects
  are to be bound together by `smbind`. 1: matrices are stacked (row
  bound), 2: matrices are column bound, 3: iterations are combined.
  Needs to be of length 1. The other dimensions are expected to match
  across the objects.

- strict:

  Logical, if consistency of the time series attributes (`"start"`,
  `"end"`, `"thin"`, and number of simulated matrices) of `simmat`
  objects are strictly enforced when binding multiple objects together
  using `smbind`. Applies only to input objects based on sequential null
  model algorithms.

- ...:

  Additional arguments supplied to algorithms. In case of `smbind` it
  can contain multiple `simmat` objects.

## Details

The purpose of the `nullmodel` function is to create an object, where
all necessary statistics of the input matrix are calculated only once.
This information is reused, but not recalculated in each step of the
simulation process done by the `simulate` method.

The `simulate` method carries out the simulation, the simulated matrices
are stored in an array. For sequential algorithms, the method updates
the state of the input `nullmodel` object. Therefore, it is possible to
do diagnostic tests on the returned `simmat` object, and make further
simulations, or use increased thinning value if desired.

The `update` method makes `burnin` steps in case of sequential
algorithms to update the status of the input model without any attempt
to return matrices. For non-sequential algorithms the method does
nothing.

`update` is the preferred way of making burnin iterations without
sampling. Alternatively, burnin can be done via the `simulate` method.
For convergence diagnostics, it is recommended to use the `simulate`
method without burnin. The input nullmodel object is updated, so further
samples can be simulated if desired without having to start the process
all over again. See Examples.

The `smbind` function can be used to combine multiple `simmat` objects.
This comes handy when null model simulations are stratified by sites
(`MARGIN = 1`) or by species (`MARGIN = 2`), or in the case when
multiple objects are returned by identical/consistent settings e.g.
during parallel computations (`MARGIN = 3`). Sanity checks are made to
ensure that combining multiple objects is sensible, but it is the user's
responsibility to check independence of the simulated matrices and the
null distribution has converged in case of sequential null model
algorithms. The `strict = FALSE` setting can relax checks regarding
start, end, and thinning values for sequential null models.

## Value

The function `nullmodel` returns an object of class `nullmodel`. It is a
set of objects sharing the same environment:

- `data`: :

  original matrix in integer mode.

- `nrow`: :

  number of rows.

- `ncol`: :

  number of columns.

- `rowSums`: :

  row sums.

- `colSums`: :

  column sums.

- `rowFreq`: :

  row frequencies (number of nonzero cells).

- `colFreq`: :

  column frequencies (number of nonzero cells).

- `totalSum`: :

  total sum.

- `fill`: :

  number of nonzero cells in the matrix.

- `commsim`: :

  the `commsim` object as a result of the `method` argument.

- `state`: :

  current state of the permutations, a matrix similar to the original.
  It is `NULL` for non-sequential algorithms.

- `iter`: :

  current number of iterations for sequential algorithms. It is `NULL`
  for non-sequential algorithms.

The `simulate` method returns an object of class `simmat`. It is an
array of simulated matrices (third dimension corresponding to `nsim`
argument).

The `update` method returns the current state (last updated matrix)
invisibly, and update the input object for sequential algorithms. For
non sequential algorithms, it returns `NULL`.

The `smbind` function returns an object of class `simmat`.

## Author

Jari Oksanen and Peter Solymos

## See also

[`commsim`](https://vegandevs.github.io/vegan/reference/commsim.md),
[`make.commsim`](https://vegandevs.github.io/vegan/reference/commsim.md),
[`permatfull`](https://vegandevs.github.io/vegan/reference/permatfull.md),
[`permatswap`](https://vegandevs.github.io/vegan/reference/permatfull.md)

## Note

Care must be taken when the input matrix only contains a single row or
column. Such input is invalid for swapping and several other methods.
This also applies to cases when the input is stratified into subsets. In
particular, subsetting can generate small or degenerate matrices that
cannot be analysed with the selected (or any) null model. These cases
are usually detected in
[`commsim`](https://vegandevs.github.io/vegan/reference/commsim.md) and
give an error. If you want to handle smoothly error cases, you should
wrap `simulate` in [`try`](https://rdrr.io/r/base/try.html) or
[`tryCatch`](https://rdrr.io/r/base/conditions.html).

## Examples

``` r
data(mite)
x <- as.matrix(mite)[1:12, 21:30]

## non-sequential nullmodel
(nm <- nullmodel(x, "r00"))
#> An object of class “nullmodel” 
#> ‘r00’ method (binary, non-sequential)
#> 12 x 10 matrix
#> 
(sm <- simulate(nm, nsim=10))
#> An object of class “simmat” 
#> ‘r00’ method (binary, non-sequential)
#> 12 x 10 matrix
#> Number of permuted matrices = 10 
#> 

## sequential nullmodel
(nm <- nullmodel(x, "swap"))
#> An object of class “nullmodel” 
#> ‘swap’ method (binary, sequential)
#> 12 x 10 matrix
#> Iterations = 0 
#> 
(sm1 <- simulate(nm, nsim=10, thin=5))
#> An object of class “simmat” 
#> ‘swap’ method (binary, sequential)
#> 12 x 10 matrix
#> Number of permuted matrices = 10 
#> Start = 5, End = 50, Thin = 5
#> 
(sm2 <- simulate(nm, nsim=10, thin=5))
#> An object of class “simmat” 
#> ‘swap’ method (binary, sequential)
#> 12 x 10 matrix
#> Number of permuted matrices = 10 
#> Start = 55, End = 100, Thin = 5
#> 

## sequential nullmodel with burnin and extra updating
(nm <- nullmodel(x, "swap"))
#> An object of class “nullmodel” 
#> ‘swap’ method (binary, sequential)
#> 12 x 10 matrix
#> Iterations = 0 
#> 
(sm1 <- simulate(nm, burnin=10, nsim=10, thin=5))
#> An object of class “simmat” 
#> ‘swap’ method (binary, sequential)
#> 12 x 10 matrix
#> Number of permuted matrices = 10 
#> Start = 15, End = 60, Thin = 5
#> 
(sm2 <- simulate(nm, nsim=10, thin=5))
#> An object of class “simmat” 
#> ‘swap’ method (binary, sequential)
#> 12 x 10 matrix
#> Number of permuted matrices = 10 
#> Start = 5, End = 50, Thin = 5
#> 

## sequential nullmodel with separate initial burnin
(nm <- nullmodel(x, "swap"))
#> An object of class “nullmodel” 
#> ‘swap’ method (binary, sequential)
#> 12 x 10 matrix
#> Iterations = 0 
#> 
nm <- update(nm, nsim=10)
(sm2 <- simulate(nm, nsim=10, thin=5))
#> An object of class “simmat” 
#> ‘swap’ method (binary, sequential)
#> 12 x 10 matrix
#> Number of permuted matrices = 10 
#> Start = 15, End = 60, Thin = 5
#> 

## combining multiple simmat objects

## stratification
nm1 <- nullmodel(x[1:6,], "r00")
sm1 <- simulate(nm1, nsim=10)
nm2 <- nullmodel(x[7:12,], "r00")
sm2 <- simulate(nm2, nsim=10)
smbind(sm1, sm2, MARGIN=1)
#> An object of class “simmat” 
#> ‘r00’ method (binary, non-sequential)
#> 12 x 10 matrix
#> Number of permuted matrices = 10 
#> 

## binding subsequent samples from sequential algorithms
## start, end, thin retained
nm <- nullmodel(x, "swap")
nm <- update(nm, nsim=10)
sm1 <- simulate(nm, nsim=10, thin=5)
sm2 <- simulate(nm, nsim=20, thin=5)
sm3 <- simulate(nm, nsim=10, thin=5)
smbind(sm3, sm2, sm1, MARGIN=3)
#> An object of class “simmat” 
#> ‘swap’ method (binary, sequential)
#> 12 x 10 matrix
#> Number of permuted matrices = 40 
#> Start = 15, End = 210, Thin = 5
#> 

## 'replicate' based usage which is similar to the output
## of 'parLapply' or 'mclapply' in the 'parallel' package
## start, end, thin are set, also noting number of chains
smfun <- function(x, burnin, nsim, thin) {
    nm <- nullmodel(x, "swap")
    nm <- update(nm, nsim=burnin)
    simulate(nm, nsim=nsim, thin=thin)
}
smlist <- replicate(3, smfun(x, burnin=50, nsim=10, thin=5), simplify=FALSE)
smbind(smlist, MARGIN=3) # Number of permuted matrices = 30
#> An object of class “simmat” 
#> ‘swap’ method (binary, sequential)
#> 12 x 10 matrix
#> Number of permuted matrices = 30 
#> Start = 55, End = 100, Thin = 5 (3 chains)
#> 

if (FALSE) { # \dontrun{
## parallel null model calculations
library(parallel)

if (.Platform$OS.type == "unix") {
## forking on Unix systems
smlist <- mclapply(1:3, function(i) smfun(x, burnin=50, nsim=10, thin=5))
smbind(smlist, MARGIN=3)
}

## socket type cluster, works on all platforms
cl <- makeCluster(3)
clusterEvalQ(cl, library(vegan))
clusterExport(cl, c("smfun", "x"))
smlist <- parLapply(cl, 1:3, function(i) smfun(x, burnin=50, nsim=10, thin=5))
stopCluster(cl)
smbind(smlist, MARGIN=3)
} # }
```

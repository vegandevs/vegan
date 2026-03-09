# Add New Points to NMDS ordination

Add new points to an existing
[`metaMDS`](https://vegandevs.github.io/vegan/reference/metaMDS.md) or
[`monoMDS`](https://vegandevs.github.io/vegan/reference/monoMDS.md)
ordination.

## Usage

``` r
MDSaddpoints(nmds, dis, neighbours = 5, maxit = 200)

dist2xy(dist, pick, type = c("xy", "xx"), invert = FALSE)
```

## Arguments

- nmds:

  Result object from
  [`metaMDS`](https://vegandevs.github.io/vegan/reference/metaMDS.md) or
  [`monoMDS`](https://vegandevs.github.io/vegan/reference/monoMDS.md).
  The configuration of points is fixed, but new points are added.

- dis:

  Rectangular non-symmetric dissimilarity matrix among new points (rows)
  and old fixed points (columns). Such matrix can be extracted from
  complete dissimilarities of both old and new points with `dist2xy`, or
  calculated with
  [`designdist2`](https://vegandevs.github.io/vegan/reference/designdist.md).

- neighbours:

  Number of nearest points used to get the starting locations for new
  points.

- maxit:

  Maximum number of iterations.

- dist:

  Input dissimilarities.

- pick:

  Indices (integers) of selected observations or a logical vector that
  is `TRUE` for picked items. The output will be in the original order
  and will not be reordered by this argument.

- type:

  `"xy"` returns rectangular data of picked against not picked
  observations, and `"xx"` a subset of symmetric dissimilarities.

- invert:

  Invert `pick`: drop elements listed.

## Value

Function return a list of class `"nmds"` (there are no other objects of
that type in vegan) with following elements

- points:

  Coordinates of added new points

- seeds:

  Starting coordinates for new points.

- deltastress:

  Change of stress with added points.

- iters:

  Number of iterations.

- cause:

  Cause of termination of iterations. Integer for convergence criteria
  in
  [`monoMDS`](https://vegandevs.github.io/vegan/reference/monoMDS.md).

## Details

Function provides an interface to `monoMDS` Fortran code to add new
points to an existing ordination that will be regarded as fixed. The
function has a similar role as `predict` functions with `newdata` in
Euclidean ordination (e.g.
[`predict.cca`](https://vegandevs.github.io/vegan/reference/predict.cca.md)).
Input data must be a rectangular matrix of distances among new added
points (rows) and all fixed old points (columns). Such matrices can be
extracted from complete dissimilarities with helper function `dist2xy`.
Function
[`designdist2`](https://vegandevs.github.io/vegan/reference/designdist.md)
can directly calculate such rectangular dissimilarity matrices between
sampling units (rows) in two matries. In addition,
[analogue](https://CRAN.R-project.org/package=analogue) has distance
function that can calculate dissimilarities among two matrices,
including functions that cannot be specified in
[`designdist2`](https://vegandevs.github.io/vegan/reference/designdist.md).

Great care is needed in preparing the dissimilarities for the input. The
dissimilarity index must be exactly the same as in the fixed ordination,
and columns must match old fixed points, and rows added new points.

## Examples

``` r
## Cross-validation: remove a point when performing NMDS and add as
## a new points
data(dune)
d <- vegdist(dune)
## remove point 3 from ordination
mod3 <- metaMDS(dist2xy(d, 3, "xx", invert = TRUE), trace=0)
## add point 3 to the result
MDSaddpoints(mod3, dist2xy(d, 3))
#> $points
#>          MDS1       MDS2
#> 3 -0.09726881 -0.4560917
#> 
#> $seed
#>            [,1]       [,2]
#> [1,] 0.01222493 -0.4353716
#> 
#> $deltastress
#> [1] 0.001811377
#> 
#> $iters
#> [1] 18
#> 
#> $cause
#> [1] 4
#> 
#> attr(,"class")
#> [1] "nmds"
## Use designdist2
d15 <- designdist(dune[1:15,])
m15 <- metaMDS(d15, trace=0)
MDSaddpoints(m15, designdist2(dune[1:15,], dune[16:20,]))
#> $points
#>          MDS1       MDS2
#> 16  0.4296366 -0.3068902
#> 17 -0.4122463  0.4118309
#> 18 -0.1135994  0.3120708
#> 19  0.1421922  0.4869128
#> 20  0.5892538  0.1751219
#> 
#> $seed
#>            [,1]        [,2]
#> [1,]  0.2832566 0.019571976
#> [2,] -0.1920994 0.137270818
#> [3,] -0.1529998 0.094374124
#> [4,] -0.0244597 0.137988191
#> [5,]  0.3212258 0.009941098
#> 
#> $deltastress
#> [1] 0.03394768
#> 
#> $iters
#> [1] 33
#> 
#> $cause
#> [1] 4
#> 
#> attr(,"class")
#> [1] "nmds"
```

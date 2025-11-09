# Connectedness of Dissimilarities

Function `distconnected` finds groups that are connected disregarding
dissimilarities that are at or above a threshold or `NA`. The function
can be used to find groups that can be ordinated together or transformed
by
[`stepacross`](https://vegandevs.github.io/vegan/reference/stepacross.md).
Function `no.shared` returns a logical dissimilarity object, where
`TRUE` means that sites have no species in common. This is a minimal
structure for `distconnected` or can be used to set missing values to
dissimilarities.

## Usage

``` r
distconnected(dis, toolong = 1, trace = TRUE)

no.shared(x)
```

## Arguments

- dis:

  Dissimilarity data inheriting from class `dist` or a an object, such
  as a matrix, that can be converted to a dissimilarity matrix.
  Functions
  [`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.md)
  and [`dist`](https://rdrr.io/r/stats/dist.html) are some functions
  producing suitable dissimilarity data.

- toolong:

  Shortest dissimilarity regarded as `NA`. The function uses a fuzz
  factor, so that dissimilarities close to the limit will be made `NA`,
  too. If `toolong = 0` (or negative), no dissimilarity is regarded as
  too long.

- trace:

  Summarize results of `distconnected`

- x:

  Community data.

## Details

Data sets are disconnected if they have sample plots or groups of sample
plots which share no species with other sites or groups of sites. Such
data sets cannot be sensibly ordinated by any unconstrained method
because these subsets cannot be related to each other. For instance,
correspondence analysis will polarize these subsets with eigenvalue 1.
Neither can such dissimilarities be transformed with
[`stepacross`](https://vegandevs.github.io/vegan/reference/stepacross.md),
because there is no path between all points, and result will contain
`NA`s. Function `distconnected` will find such subsets in dissimilarity
matrices. The function will return a grouping vector that can be used
for sub-setting the data. If data are connected, the result vector will
be all \\1\\s. The connectedness between two points can be defined
either by a threshold `toolong` or using input dissimilarities with
`NA`s.

Function `no.shared` returns a `dist` structure having value `TRUE` when
two sites have nothing in common, and value `FALSE` when they have at
least one shared species. This is a minimal structure that can be
analysed with `distconnected`. The function can be used to select
dissimilarities with no shared species in indices which do not have a
fixed upper limit.

Function `distconnected` uses depth-first search (Sedgewick 1990).

## Value

Function `distconnected` returns a vector for observations using
integers to identify connected groups. If the data are connected, values
will be all `1`. Function `no.shared` returns an object of class
[`dist`](https://rdrr.io/r/stats/dist.html).

## References

Sedgewick, R. (1990). *Algorithms in C*. Addison Wesley.

## Author

Jari Oksanen

## See also

[`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.md) or
[`dist`](https://rdrr.io/r/stats/dist.html) for getting dissimilarities,
[`stepacross`](https://vegandevs.github.io/vegan/reference/stepacross.md)
for a case where you may need `distconnected`, and for connecting points
[`spantree`](https://vegandevs.github.io/vegan/reference/spantree.md).

## Examples

``` r
## There are no disconnected data in vegan, and the following uses an
## extremely low threshold limit for connectedness. This is for
## illustration only, and not a recommended practice.
data(dune)
dis <- vegdist(dune)
gr <- distconnected(dis, toolong=0.4)
#> Connectivity of distance matrix with threshold dissimilarity 0.4 
#> Data are disconnected: 6 groups
#> Groups sizes
#>  1  2  3  4  5  6 
#>  1 11  2  4  1  1 
# Make sites with no shared species as NA in Manhattan dissimilarities
dis <- vegdist(dune, "manhattan")
is.na(dis) <- no.shared(dune)
```

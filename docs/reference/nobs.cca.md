# Extract the Number of Observations from a vegan Fit.

Extract the number of ‘observations’ from a vegan model fit.

## Usage

``` r
# S3 method for class 'cca'
nobs(object, ...)
```

## Arguments

- object:

  A fitted model object.

- ...:

  Further arguments to be passed to methods.

## Details

Function `nobs` is generic in R, and vegan provides methods for objects
from
[`betadisper`](https://vegandevs.github.io/vegan/reference/betadisper.md),
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md) and other
related methods,
[`CCorA`](https://vegandevs.github.io/vegan/reference/CCorA.md),
[`decorana`](https://vegandevs.github.io/vegan/reference/decorana.md),
[`isomap`](https://vegandevs.github.io/vegan/reference/isomap.md),
[`metaMDS`](https://vegandevs.github.io/vegan/reference/metaMDS.md),
[`pcnm`](https://vegandevs.github.io/vegan/reference/pcnm.md),
[`procrustes`](https://vegandevs.github.io/vegan/reference/procrustes.md),
[`radfit`](https://vegandevs.github.io/vegan/reference/radfit.md),
[`varpart`](https://vegandevs.github.io/vegan/reference/varpart.md) and
[`wcmdscale`](https://vegandevs.github.io/vegan/reference/wcmdscale.md).

## Value

A single number, normally an integer, giving the number of observations.

## Author

Jari Oksanen

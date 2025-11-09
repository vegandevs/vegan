# Vegetation and environment in lichen pastures

The `varespec` data frame has 24 rows and 44 columns. Columns are
estimated cover values of 44 species. The variable names are formed from
the scientific names, and are self explanatory for anybody familiar with
the vegetation type. The `varechem` data frame has 24 rows and 14
columns, giving the soil characteristics of the very same sites as in
the `varespec` data frame. The chemical measurements have obvious names.
`Baresoil` gives the estimated cover of bare soil, `Humdepth` the
thickness of the humus layer.

## Usage

``` r
data(varechem)
       data(varespec)
```

## References

Väre, H., Ohtonen, R. and Oksanen, J. (1995) Effects of reindeer grazing
on understorey vegetation in dry Pinus sylvestris forests. *Journal of
Vegetation Science* 6, 523–530.

## Examples

``` r
data(varespec)
data(varechem)
```

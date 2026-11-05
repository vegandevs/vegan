# Norwegian Continental Shelf Macrofauna Data

Counts of 809 benthic soft-sediment macrofaunal species at 101 sites
from five areas along the Norwegian continental shelf. Anderson et al.
(2017) used presence/absence versions of these data to illustrate the
modified PERMANOVA statistic for heterogeneous dispersions.

## Usage

``` r
data("norwegianShelf")
```

## Format

A list with three components:

- comm:

  An integer matrix with 101 sites as rows and 809 species as columns.

- area:

  A factor assigning the 101 sites to five sampling areas.

- site:

  Integer site identifiers from 1 to 101.

## Details

The published presence/absence analysis is obtained with
`1L * (norwegianShelf$comm > 0L)`.

## Source

Supporting Data S1 for Anderson et al. (2017),
[doi:10.1111/anzs.12176](https://doi.org/10.1111/anzs.12176) .

## References

Anderson, M.J., Walsh, D.C.I., Clarke, K.R., Gorley, R.N. and Guerra-
Castro, E. (2017). Some solutions to the multivariate Behrens–Fisher
problem for dissimilarity-based analyses. *Australian & New Zealand
Journal of Statistics* 59, 57–79.

## Examples

``` r
data("norwegianShelf")
str(norwegianShelf, max.level = 1)
#> List of 3
#>  $ comm: int [1:101, 1:809] 1 0 0 1 0 0 0 0 0 0 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ area: Factor w/ 5 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ site: int [1:101] 1 2 3 4 5 6 7 8 9 10 ...
```

# Morisita index of intraspecific aggregation

Calculates the Morisita index of dispersion, standardized index values,
and the so called clumpedness and uniform indices.

## Usage

``` r
dispindmorisita(x, unique.rm = FALSE, crit = 0.05, na.rm = FALSE)
```

## Arguments

- x:

  community data matrix, with sites (samples) as rows and species as
  columns.

- unique.rm:

  logical, if `TRUE`, unique species (occurring in only one sample) are
  removed from the result.

- crit:

  two-sided p-value used to calculate critical Chi-squared values.

- na.rm:

  logical. Should missing values (including `NaN`) be omitted from the
  calculations?

## Details

The Morisita index of dispersion is defined as (Morisita 1959, 1962):

`Imor = n * (sum(xi^2) - sum(xi)) / (sum(xi)^2 - sum(xi))`

where \\xi\\ is the count of individuals in sample \\i\\, and \\n\\ is
the number of samples (\\i = 1, 2, \ldots, n\\). \\Imor\\ has values
from 0 to \\n\\. In uniform (hyperdispersed) patterns its value falls
between 0 and 1, in clumped patterns it falls between 1 and \\n\\. For
increasing sample sizes (i.e. joining neighbouring quadrats), \\Imor\\
goes to \\n\\ as the quadrat size approaches clump size. For random
patterns, \\Imor = 1\\ and counts in the samples follow Poisson
frequency distribution.

The deviation from random expectation (null hypothesis) can be tested
using critical values of the Chi-squared distribution with \\n-1\\
degrees of freedom. Confidence intervals around 1 can be calculated by
the clumped \\Mclu\\ and uniform \\Muni\\ indices (Hairston et al. 1971,
Krebs 1999) (Chi2Lower and Chi2Upper refers to e.g. 0.025 and 0.975
quantile values of the Chi-squared distribution with \\n-1\\ degrees of
freedom, respectively, for `crit = 0.05`):

`Mclu = (Chi2Lower - n + sum(xi)) / (sum(xi) - 1)`

`Muni = (Chi2Upper - n + sum(xi)) / (sum(xi) - 1)`

Smith-Gill (1975) proposed scaling of Morisita index from \[0, n\]
interval into \[-1, 1\], and setting up -0.5 and 0.5 values as
confidence limits around random distribution with rescaled value 0. To
rescale the Morisita index, one of the following four equations apply to
calculate the standardized index \\Imst\\:

\(a\) `Imor >= Mclu > 1`: `Imst = 0.5 + 0.5 (Imor - Mclu) / (n - Mclu)`,

\(b\) `Mclu > Imor >= 1`: `Imst = 0.5 (Imor - 1) / (Mclu - 1)`,

\(c\) `1 > Imor > Muni`: `Imst = -0.5 (Imor - 1) / (Muni - 1)`,

\(d\) `1 > Muni > Imor`: `Imst = -0.5 + 0.5 (Imor - Muni) / Muni`.

## Value

Returns a data frame with as many rows as the number of columns in the
input data, and with four columns. Columns are: `imor` the
unstandardized Morisita index, `mclu` the clumpedness index, `muni` the
uniform index, `imst` the standardized Morisita index, `pchisq` the
Chi-squared based probability for the null hypothesis of random
expectation.

## References

Morisita, M. 1959. Measuring of the dispersion of individuals and
analysis of the distributional patterns. *Mem. Fac. Sci. Kyushu Univ.
Ser. E* 2, 215–235.

Morisita, M. 1962. Id-index, a measure of dispersion of individuals.
*Res. Popul. Ecol.* 4, 1–7.

Smith-Gill, S. J. 1975. Cytophysiological basis of disruptive pigmentary
patterns in the leopard frog, *Rana pipiens*. II. Wild type and mutant
cell specific patterns. *J. Morphol.* 146, 35–54.

Hairston, N. G., Hill, R. and Ritte, U. 1971. The interpretation of
aggregation patterns. In: Patil, G. P., Pileou, E. C. and Waters, W. E.
eds. *Statistical Ecology 1: Spatial Patterns and Statistical
Distributions*. Penn. State Univ. Press, University Park.

Krebs, C. J. 1999. *Ecological Methodology*. 2nd ed. Benjamin Cummings
Publishers.

## Author

Péter Sólymos, <solymos@ualberta.ca>

## Note

A common error found in several papers is that when standardizing as in
the case (b), the denominator is given as `Muni - 1`. This results in a
hiatus in the \[0, 0.5\] interval of the standardized index. The root of
this typo is the book of Krebs (1999), see the Errata for the book (Page
217, currently
https://www.zoology.ubc.ca/~krebs/downloads/errors_2nd_printing.pdf).

## Examples

``` r
data(dune)
x <- dispindmorisita(dune)
x
#>                imor      mclu         muni        imst       pchisq
#> Achimill  2.1666667  1.923488  0.327101099  0.50672636 9.157890e-03
#> Agrostol  1.8085106  1.294730  0.785245032  0.51373357 1.142619e-05
#> Airaprae  8.0000000  4.463082 -1.523370880  0.61382303 3.571702e-04
#> Alopgeni  2.5396825  1.395781  0.711614757  0.53074307 3.024441e-08
#> Anthodor  2.6666667  1.692616  0.495325824  0.52660266 5.897217e-05
#> Bellpere  2.0512821  2.154361  0.158876373  0.45535255 3.451547e-02
#> Bromhord  3.2380952  1.989452  0.279036892  0.53466422 1.170437e-04
#> Chenalbu        NaN       Inf         -Inf         NaN          NaN
#> Cirsarve 20.0000000 14.852327 -9.093483518  1.00000000 5.934709e-03
#> Comapalu  6.6666667  5.617442 -2.364494506  0.53647558 1.055552e-02
#> Eleopalu  3.7333333  1.577180  0.579438187  0.55851854 2.958285e-10
#> Elymrepe  2.7692308  1.554093  0.596260659  0.53293787 1.180195e-06
#> Empenigr 20.0000000 14.852327 -9.093483518  1.00000000 5.934709e-03
#> Hyporadi  6.6666667  2.731541 -0.261685440  0.61393969 7.832274e-07
#> Juncarti  3.1372549  1.814843  0.406265675  0.53635966 2.066336e-05
#> Juncbufo  4.1025641  2.154361  0.158876373  0.55458486 1.503205e-05
#> Lolipere  1.5849970  1.243023  0.822921342  0.50911591 5.873839e-05
#> Planlanc  2.4615385  1.554093  0.596260659  0.52459747 1.921730e-05
#> Poaprat   1.1702128  1.294730  0.785245032  0.28876015 1.046531e-01
#> Poatriv   1.4644137  1.223425  0.837201879  0.50641728 2.747301e-04
#> Ranuflam  2.4175824  2.065564  0.223578191  0.50981405 7.010483e-03
#> Rumeacet  3.9215686  1.814843  0.406265675  0.55792432 1.530085e-07
#> Sagiproc  2.4210526  1.729070  0.468764025  0.51893672 4.956394e-04
#> Salirepe  5.8181818  2.385233 -0.009348352  0.59744520 2.687397e-07
#> Scorautu  0.9643606  1.261365  0.809556915 -0.09356972 5.823404e-01
#> Trifprat  6.6666667  2.731541 -0.261685440  0.61393969 7.832274e-07
#> Trifrepe  1.2210916  1.301138  0.780576445  0.36709402 6.335449e-02
#> Vicilath  3.3333333  5.617442 -2.364494506  0.25266513 1.301890e-01
#> Bracruta  1.1904762  1.288590  0.789719093  0.33001160 8.071762e-02
#> Callcusp  5.3333333  2.539147 -0.121498169  0.58001287 7.982634e-06
y <- dispindmorisita(dune, unique.rm = TRUE)
y
#>               imor     mclu         muni        imst       pchisq
#> Achimill 2.1666667 1.923488  0.327101099  0.50672636 9.157890e-03
#> Agrostol 1.8085106 1.294730  0.785245032  0.51373357 1.142619e-05
#> Airaprae 8.0000000 4.463082 -1.523370880  0.61382303 3.571702e-04
#> Alopgeni 2.5396825 1.395781  0.711614757  0.53074307 3.024441e-08
#> Anthodor 2.6666667 1.692616  0.495325824  0.52660266 5.897217e-05
#> Bellpere 2.0512821 2.154361  0.158876373  0.45535255 3.451547e-02
#> Bromhord 3.2380952 1.989452  0.279036892  0.53466422 1.170437e-04
#> Comapalu 6.6666667 5.617442 -2.364494506  0.53647558 1.055552e-02
#> Eleopalu 3.7333333 1.577180  0.579438187  0.55851854 2.958285e-10
#> Elymrepe 2.7692308 1.554093  0.596260659  0.53293787 1.180195e-06
#> Hyporadi 6.6666667 2.731541 -0.261685440  0.61393969 7.832274e-07
#> Juncarti 3.1372549 1.814843  0.406265675  0.53635966 2.066336e-05
#> Juncbufo 4.1025641 2.154361  0.158876373  0.55458486 1.503205e-05
#> Lolipere 1.5849970 1.243023  0.822921342  0.50911591 5.873839e-05
#> Planlanc 2.4615385 1.554093  0.596260659  0.52459747 1.921730e-05
#> Poaprat  1.1702128 1.294730  0.785245032  0.28876015 1.046531e-01
#> Poatriv  1.4644137 1.223425  0.837201879  0.50641728 2.747301e-04
#> Ranuflam 2.4175824 2.065564  0.223578191  0.50981405 7.010483e-03
#> Rumeacet 3.9215686 1.814843  0.406265675  0.55792432 1.530085e-07
#> Sagiproc 2.4210526 1.729070  0.468764025  0.51893672 4.956394e-04
#> Salirepe 5.8181818 2.385233 -0.009348352  0.59744520 2.687397e-07
#> Scorautu 0.9643606 1.261365  0.809556915 -0.09356972 5.823404e-01
#> Trifprat 6.6666667 2.731541 -0.261685440  0.61393969 7.832274e-07
#> Trifrepe 1.2210916 1.301138  0.780576445  0.36709402 6.335449e-02
#> Vicilath 3.3333333 5.617442 -2.364494506  0.25266513 1.301890e-01
#> Bracruta 1.1904762 1.288590  0.789719093  0.33001160 8.071762e-02
#> Callcusp 5.3333333 2.539147 -0.121498169  0.58001287 7.982634e-06
dim(x) ## with unique species
#> [1] 30  5
dim(y) ## unique species removed
#> [1] 27  5
```

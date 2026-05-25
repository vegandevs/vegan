# Species tolerances and sample heterogeneities

Species tolerances and sample heterogeneities.

## Usage

``` r
tolerance(x, ...)

# S3 method for class 'cca'
tolerance(x, choices = 1:2, which = c("species","sites"),
          scaling = "species", useN2 = TRUE, hill = FALSE, ...)

# S3 method for class 'decorana'
tolerance(x, data, choices = 1:4,
          which = c("sites", "species"), useN2 = TRUE, ...)
```

## Details

Function to compute species tolerances and site heterogeneity measures
from unimodal ordinations (CCA & CA). Implements Eq 6.47 and 6.48 from
the Canoco 4.5 Reference Manual (pages 178–179).

Function `wascores` with `stdev = TRUE` uses the same algebra, but bases
the standard deviations on weighted averages scores instead of linear
combinations scores of `tolerance`.

## Value

Matrix of tolerances/heterogeneities with some additional attributes:
`which`, `scaling`, and `N2`, the latter of which will be `NA` if
`useN2 = FALSE` or `N2` could not be estimated.

## Author

Gavin L. Simpson and Jari Oksanen (`decorana` method).

## Arguments

- x:

  object of class `"cca"`.

- choices:

  numeric; which ordination axes to compute tolerances and
  heterogeneities for. Defaults to axes 1 and 2.

- which:

  character; one of `"species"` or `"sites"`, indicating whether species
  tolerances or sample heterogeneities respectively are computed.

- scaling:

  character or numeric; the ordination scaling to use. See
  [`scores.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md)
  for details.

- hill:

  logical; if `scaling` is a character, these control whether Hill's
  scaling is used for (C)CA respectively. See
  [`scores.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md)
  for details.

- useN2:

  logical; should the bias in the tolerances / heterogeneities be
  reduced via scaling by Hill's N2?

- data:

  Original input data used in
  [`decorana`](https://vegandevs.github.io/vegan/reference/decorana.md).
  If missing, the function tries to get the same data as used in
  `decorana` call.

- ...:

  arguments passed to other methods.

## Examples

``` r
data(dune)
data(dune.env)
mod <- cca(dune ~ ., data = dune.env)
#> 
#> Some constraints or conditions were aliased because they were redundant. This
#> can happen if terms are constant or linearly dependent (collinear): ‘Manure^4’

## defaults to species tolerances
tolerance(mod)
#> 
#> Species Tolerance
#> 
#> Scaling: 2
#> 
#>                CCA1      CCA2
#> Achimill 0.32968099 0.9241988
#> Agrostol 0.93670069 0.9238455
#> Airaprae 1.04694096 0.5889849
#> Alopgeni 0.72227472 0.3760138
#> Anthodor 1.00596787 0.8338212
#> Bellpere 0.32891011 0.9962790
#> Bromhord 0.27740999 0.6236199
#> Chenalbu 0.00000000 0.0000000
#> Cirsarve 0.00000000 0.0000000
#> Comapalu 0.47185632 0.8029414
#> Eleopalu 0.50344134 0.9384960
#> Elymrepe 0.35119963 0.5642491
#> Empenigr 0.00000000 0.0000000
#> Hyporadi 1.05840696 0.7523003
#> Juncarti 0.78397702 1.0686743
#> Juncbufo 0.69275956 0.6180830
#> Lolipere 0.51006235 0.8278177
#> Planlanc 0.36040676 0.6962294
#> Poaprat  0.58184277 0.9547104
#> Poatriv  0.78695928 0.7433503
#> Ranuflam 0.56576326 1.1725628
#> Rumeacet 0.58715663 0.8751491
#> Sagiproc 0.70922180 1.1153129
#> Salirepe 0.98530179 0.1077917
#> Scorautu 1.04355761 1.0724439
#> Trifprat 0.03045846 0.3651949
#> Trifrepe 1.21543364 0.9115613
#> Vicilath 0.24853962 0.6194084
#> Bracruta 1.03787313 1.0958331
#> Callcusp 0.57882025 1.0418623
#> 

## sample heterogeneities for CCA axes 1:6
tolerance(mod, which = "sites", choices = 1:6)
#> 
#> Sample Heterogeneity
#> 
#> Scaling: 2
#> 
#>         CCA1      CCA2      CCA3      CCA4      CCA5      CCA6
#> 1  0.2350112 0.8611530 1.7964571 0.4445499 2.4235732 0.5496289
#> 2  0.7100754 0.4136311 0.8151643 0.6311751 1.0467901 0.2514646
#> 3  0.5076492 0.7279717 0.8306874 0.5590739 0.3904998 0.9162012
#> 4  0.5955037 0.6901907 0.7931255 0.4873638 0.3966068 0.8700581
#> 5  0.6001048 0.5614830 1.1481560 0.3569604 0.4423909 1.9420043
#> 6  0.7272637 0.6867342 1.6068628 0.7778498 0.9187843 0.4938865
#> 7  0.6478967 0.4993262 0.7207318 0.3817131 0.4130713 0.7228173
#> 8  0.8563491 0.5498552 0.4217718 0.3370226 0.3013276 0.9535190
#> 9  0.5599722 0.7399384 0.4170304 1.0535541 1.4612437 0.7626183
#> 10 0.5210280 0.5806978 0.5856634 0.4174860 1.8559344 0.8890262
#> 11 0.4489323 0.6016877 0.3317371 1.8780211 1.2965939 2.1953737
#> 12 0.4948094 1.1084494 0.5226746 1.5064446 0.5703077 1.1561020
#> 13 0.6998985 0.8859365 0.4215474 0.8582272 0.5673698 0.5186678
#> 14 1.5925779 0.6747926 0.8927360 1.6798300 0.3480218 0.1575892
#> 15 1.0107648 0.5294221 1.0975629 1.7632888 0.2240900 0.3727240
#> 16 0.8031479 0.6058313 0.4871527 0.4227451 0.5341256 0.6990815
#> 17 0.5936276 1.5142792 0.5137979 1.0224938 1.7931775 0.6261853
#> 18 0.5689409 1.4067575 0.6398557 0.4983399 0.4364791 0.6590394
#> 19 1.1330387 0.9816332 1.1242398 0.7238920 0.5577662 0.7036044
#> 20 0.6737757 1.4458326 1.4380928 1.0959027 0.4142423 0.5332460
#> 
## average should be 1 with scaling = "sites", hill = TRUE
tol <- tolerance(mod, which = "sites", scaling = "sites", hill = TRUE,
   choices = 1:4)
colMeans(tol)
#>     CCA1     CCA2     CCA3     CCA4 
#> 1.059199 1.048823 1.000551 1.077612 
apply(tol, 2, sd)
#>      CCA1      CCA2      CCA3      CCA4 
#> 0.3174462 0.2793521 0.3714540 0.2681931 
## Rescaling tries to set all tolerances to 1
tol <- tolerance(decorana(dune))
colMeans(tol)
#>      DCA1      DCA2      DCA3      DCA4 
#> 0.9817657 0.9249544 0.9444811 0.9821666 
apply(tol, 2, sd)
#>      DCA1      DCA2      DCA3      DCA4 
#> 0.1977777 0.3204058 0.2646872 0.1210543 
```

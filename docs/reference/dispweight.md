# Dispersion-based weighting of species counts

Transform abundance data downweighting species that are overdispersed to
the Poisson error.

## Usage

``` r
dispweight(comm, groups, nsimul = 999, nullmodel = "c0_ind",
    plimit = 0.05)
gdispweight(formula, data, plimit = 0.05)
# S3 method for class 'dispweight'
summary(object, ...)
```

## Arguments

- comm:

  Community data matrix.

- groups:

  Factor describing the group structure. If missing, all sites are
  regarded as belonging to one group. `NA` values are not allowed.

- nsimul:

  Number of simulations.

- nullmodel:

  The
  [`nullmodel`](https://vegandevs.github.io/vegan/reference/nullmodel.md)
  used in
  [`commsim`](https://vegandevs.github.io/vegan/reference/commsim.md)
  within `groups`. The default follows Clarke et al. (2006).

- plimit:

  Downweight species if their \\p\\-value is at or below this limit.

- formula, data:

  Formula where the left-hand side is the community data frame and
  right-hand side gives the explanatory variables. The explanatory
  variables are found in the data frame given in `data` or in the parent
  frame.

- object:

  Result object from `dispweight` or `gdispweight`.

- ...:

  Other parameters passed to functions.

## Details

The dispersion index (\\D\\) is calculated as ratio between variance and
expected value for each species. If the species abundances follow
Poisson distribution, expected dispersion is \\E(D) = 1\\, and if \\D \>
1\\, the species is overdispersed. The inverse \\1/D\\ can be used to
downweight species abundances. Species are only downweighted when
overdispersion is judged to be statistically significant (Clarke et al.
2006).

Function `dispweight` implements the original procedure of Clarke et al.
(2006). Only one factor can be used to group the sites and to find the
species means. The significance of overdispersion is assessed freely
distributing individuals of each species within factor levels. This is
achieved by using
[`nullmodel`](https://vegandevs.github.io/vegan/reference/nullmodel.md)
`"c0_ind"` (which accords to Clarke et al. 2006), but other nullmodels
can be used, though they may not be meaningful (see
[`commsim`](https://vegandevs.github.io/vegan/reference/commsim.md) for
alternatives). If a species is absent in some factor level, the whole
level is ignored in calculation of overdispersion, and the number of
degrees of freedom can vary among species. The reduced number of degrees
of freedom is used as a divisor for overdispersion \\D\\, and such
species have higher dispersion and hence lower weights in
transformation.

Function `gdispweight` is a generalized parametric version of
`dispweight`. The function is based on
[`glm`](https://rdrr.io/r/stats/glm.html) with
[`quasipoisson`](https://rdrr.io/r/stats/family.html) error
[`family`](https://rdrr.io/r/stats/family.html). Any
[`glm`](https://rdrr.io/r/stats/glm.html) model can be used, including
several factors or continuous covariates. Function `gdispweight` uses
the same test statistic as `dispweight` (Pearson Chi-square), but it
does not ignore factor levels where species is absent, and the number of
degrees of freedom is equal for all species. Therefore transformation
weights can be higher than in `dispweight`. The `gdispweight` function
evaluates the significance of overdispersion parametrically from
Chi-square distribution
([`pchisq`](https://rdrr.io/r/stats/Chisquare.html)).

Functions `dispweight` and `gdispweight` transform data, but they add
information on overdispersion and weights as attributes of the result.
The `summary` can be used to extract and print that information.

## Value

Function returns transformed data with the following new attributes:

- D:

  Dispersion statistic.

- df:

  Degrees of freedom for each species.

- p:

  \\p\\-value of the Dispersion statistic \\D\\.

- weights:

  weights applied to community data.

- nsimul:

  Number of simulations used to assess the \\p\\-value, or `NA` when
  simulations were not performed.

- nullmodel:

  The name of
  [`commsim`](https://vegandevs.github.io/vegan/reference/commsim.md)
  null model, or `NA` when simulations were not performed.

## References

Clarke, K. R., M. G. Chapman, P. J. Somerfield, and H. R. Needham. 2006.
Dispersion-based weighting of species counts in assemblage analyses.
*Marine Ecology Progress Series*, 320, 11–27.

## Author

Eduard Szöcs <eduardszoesc@gmail.com> wrote the original `dispweight`,
Jari Oksanen significantly modified the code, provided support functions
and developed `gdispweight`.

## Examples

``` r
data(mite, mite.env)
## dispweight and its summary
mite.dw <- with(mite.env, dispweight(mite, Shrub, nsimul = 99))
## IGNORE_RDIFF_BEGIN
summary(mite.dw)
#>          Dispersion    Weight Df Pr(Disp.)   
#> Brachy       9.6908 0.1031909 67      0.01 **
#> PHTH         3.2809 0.3047900 49      0.01 **
#> HPAV         6.5263 0.1532264 67      0.01 **
#> RARD         6.0477 0.1653525 49      0.01 **
#> SSTR         2.2619 0.4421053 49      0.01 **
#> Protopl      5.4229 0.1844031 49      0.01 **
#> MEGR         4.5354 0.2204860 67      0.01 **
#> MPRO         1.2687 0.7882353 67      0.05 * 
#> TVIE         2.5956 0.3852706 67      0.01 **
#> HMIN        10.0714 0.0992906 67      0.01 **
#> HMIN2        7.5674 0.1321466 49      0.01 **
#> NPRA         2.6743 0.3739344 67      0.01 **
#> TVEL         9.6295 0.1038474 49      0.01 **
#> ONOV        11.3628 0.0880064 67      0.01 **
#> SUCT         8.7372 0.1144533 67      0.01 **
#> LCIL       129.4436 0.0077254 67      0.01 **
#> Oribatl1     4.1250 0.2424248 67      0.01 **
#> Ceratoz1     1.7150 0.5830768 67      0.01 **
#> PWIL         2.2943 0.4358538 67      0.01 **
#> Galumna1     2.8777 0.3474943 49      0.01 **
#> Stgncrs2     3.8242 0.2614953 49      0.01 **
#> HRUF         1.7575 0.5690021 67      0.01 **
#> Trhypch1    14.9225 0.0670128 67      0.01 **
#> PPEL         1.3628 1.0000000 49      0.08 . 
#> NCOR         2.5875 0.3864771 67      0.01 **
#> SLAT         2.7857 0.3589744 49      0.01 **
#> FSET         4.8901 0.2044944 49      0.01 **
#> Lepidzts     1.6577 0.6032360 49      0.01 **
#> Eupelops     1.4611 0.6844033 67      0.02 * 
#> Miniglmn     1.6505 0.6058733 49      0.01 **
#> LRUG        12.0658 0.0828790 67      0.01 **
#> PLAG2        3.2403 0.3086090 67      0.01 **
#> Ceratoz3     3.5125 0.2846947 67      0.01 **
#> Oppiminu     3.1680 0.3156525 67      0.01 **
#> Trimalc2    10.5927 0.0944046 67      0.01 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Based on 99 simulations on 'c0_ind' nullmodel
## IGNORE_RDIFF_END
## generalized dispersion weighting
mite.dw <- gdispweight(mite ~ Shrub + WatrCont, data = mite.env)
rda(mite.dw ~ Shrub + WatrCont, data = mite.env)
#> 
#> Call: rda(formula = mite.dw ~ Shrub + WatrCont, data = mite.env)
#> 
#>               Inertia Proportion Rank
#> Total         38.1640     1.0000     
#> Constrained    9.2129     0.2414    3
#> Unconstrained 28.9511     0.7586   35
#> 
#> Inertia is variance
#> 
#> Eigenvalues for constrained axes:
#>  RDA1  RDA2  RDA3 
#> 7.986 0.748 0.480 
#> 
#> Eigenvalues for unconstrained axes:
#>   PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8 
#> 5.886 3.634 2.791 2.592 1.932 1.573 1.210 1.078 
#> (Showing 8 of 35 unconstrained eigenvalues)
#> 
```

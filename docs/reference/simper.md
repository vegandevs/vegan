# Similarity Percentages

Discriminating species between two groups using Bray-Curtis
dissimilarities

## Usage

``` r
simper(comm, group, permutations = 999, parallel = 1, ...)
# S3 method for class 'simper'
summary(object, ordered = TRUE,
    digits = max(3,getOption("digits") - 3), ...)
```

## Arguments

- comm:

  Community data.

- group:

  Factor describing the group structure. If this is missing or has only
  one level, contributions are estimated for non-grouped data and
  dissimilarities only show the overall heterogeneity in species
  abundances.

- permutations:

  a list of control values for the permutations as returned by the
  function [`how`](https://rdrr.io/pkg/permute/man/how.html), or the
  number of permutations required, or a permutation matrix where each
  row gives the permuted indices.

- object:

  an object returned by `simper`.

- ordered:

  Logical; Should the species be ordered by their average contribution?

- digits:

  Number of digits in output.

- parallel:

  Number of parallel processes or a predefined socket cluster. With
  `parallel = 1` uses ordinary, non-parallel processing. (Not yet
  implemented).

- ...:

  Parameters passed to other functions. In `simper` the extra parameters
  are passed to
  [`shuffleSet`](https://rdrr.io/pkg/permute/man/shuffleSet.html) if
  permutations are used.

## Details

Similarity percentage, `simper` (Clarke 1993) is based on the
decomposition of Bray-Curtis dissimilarity index (see
[`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.md),
[`designdist`](https://vegandevs.github.io/vegan/reference/designdist.md)).
The contribution of individual species \\i\\ to the overall Bray-Curtis
dissimilarity \\d\_{jk}\\ is given by

\$\$d\_{ijk} = \frac{\|x\_{ij}-x\_{ik}\|}{\sum\_{i=1}^S
(x\_{ij}+x\_{ik})}\$\$

where \\x\\ is the abundance of species \\i\\ in sampling units \\j\\
and \\k\\. The overall index is the sum of the individual contributions
over all \\S\\ species \\d\_{jk}=\sum\_{i=1}^S d\_{ijk}\\.

The `simper` functions performs pairwise comparisons of groups of
sampling units and finds the contribution of each species to the average
between-group Bray-Curtis dissimilarity. Although the method is called
“Similarity Percentages”, it really studied dissimilarities instead of
similarities (Clarke 1993).

The function displays most important species for each pair of `groups`.
These species contribute at least to 70 % of the differences between
groups. The function returns much more extensive results (including all
species) which can be accessed directly from the result object (see
section Value). Function `summary` transforms the result to a list of
data frames. With argument `ordered = TRUE` the data frames also include
the cumulative contributions and are ordered by species contribution.

The results of `simper` can be very difficult to interpret and they are
often misunderstood even in publications. The method gives the
contribution of each species to overall dissimilarities, but these are
caused by variation in species abundances, and only partly by
differences among groups. Even if you make groups that are copies of
each other, the method will single out species with high contribution,
but these are not contributions to non-existing between-group
differences but to random noise variation in species abundances. The
most abundant species usually have highest variances, and they have high
contributions even when they do not differ among groups. Permutation
tests study the differences among groups, and they can be used to find
out the species for which the differences among groups is an important
component of their contribution to dissimilarities. Analysis without
`group` argument will find species contributions to the average overall
dissimilarity among sampling units. These non-grouped contributions can
be compared to grouped contributions to see how much added value the
grouping has for each species.

## Value

A list of class `"simper"` with following items:

- species:

  The species names.

- average:

  Species contribution to average between-group dissimilarity.

- overall:

  The average between-group dissimilarity. This is the sum of the item
  `average`.

- sd:

  Standard deviation of contribution.

- ratio:

  Average to sd ratio.

- ava, avb:

  Average abundances per group.

- ord:

  An index vector to order vectors by their contribution or order
  `cusum` back to the original data order.

- cusum:

  Ordered cumulative contribution. These are based on item `average`,
  but they sum up to total 1.

- p:

  Permutation \\p\\-value. Probability of getting a larger or equal
  average contribution in random permutation of the group factor. These
  area only available if `permutations` were used (default: not
  calculated).

## See also

Function
[`meandist`](https://vegandevs.github.io/vegan/reference/mrpp.md) shows
the average between-group dissimilarities (as well as the within-group
dissimilarities).

## Author

Eduard Szöcs and Jari Oksanen.

## References

Clarke, K.R. 1993. Non-parametric multivariate analyses of changes in
community structure. *Australian Journal of Ecology*, 18, 117–143.

## Examples

``` r
data(dune)
data(dune.env)
(sim <- with(dune.env, simper(dune, Management, permutations = 99)))
#> cumulative contributions of most influential species:
#> 
#> $SF_BF
#>   Agrostol   Alopgeni   Lolipere   Trifrepe    Poatriv   Scorautu   Bromhord 
#> 0.09824271 0.18254830 0.25956958 0.33367870 0.40734444 0.47729205 0.53120026 
#>   Achimill   Planlanc   Elymrepe   Bracruta 
#> 0.57946526 0.62522255 0.67016196 0.71098133 
#> 
#> $SF_HF
#>   Agrostol   Alopgeni   Lolipere   Planlanc   Rumeacet   Elymrepe    Poatriv 
#> 0.08350879 0.16534834 0.23934930 0.30843624 0.37716139 0.43334492 0.48351753 
#>   Bracruta   Eleopalu    Poaprat   Anthodor   Sagiproc   Trifprat 
#> 0.52804045 0.57205850 0.61423981 0.65549838 0.69628951 0.73696831 
#> 
#> $SF_NM
#>   Poatriv  Alopgeni  Agrostol  Lolipere  Eleopalu   Poaprat  Bracruta  Elymrepe 
#> 0.1013601 0.1935731 0.2667383 0.3377578 0.3999419 0.4526707 0.5044725 0.5505643 
#>  Scorautu  Trifrepe  Sagiproc  Salirepe 
#> 0.5926117 0.6320111 0.6712478 0.7091528 
#> 
#> $BF_HF
#>   Rumeacet    Poatriv   Planlanc   Bromhord   Lolipere   Elymrepe   Trifrepe 
#> 0.08163219 0.15193797 0.21918333 0.27967181 0.33969561 0.39843338 0.45298204 
#>   Anthodor   Achimill   Bracruta   Alopgeni   Trifprat   Juncarti 
#> 0.50276849 0.55222648 0.60021994 0.64584333 0.69126471 0.73366621 
#> 
#> $BF_NM
#>  Lolipere   Poatriv   Poaprat  Trifrepe  Bromhord  Bracruta  Eleopalu  Agrostol 
#> 0.1242718 0.1992126 0.2711756 0.3414609 0.3958520 0.4448077 0.4910724 0.5369083 
#>  Achimill  Scorautu  Anthodor  Planlanc 
#> 0.5823926 0.6253645 0.6638182 0.7012577 
#> 
#> $HF_NM
#>    Poatriv   Lolipere   Rumeacet    Poaprat   Planlanc   Bracruta   Eleopalu 
#> 0.09913221 0.17468460 0.23917190 0.29701331 0.35469313 0.40365488 0.44804851 
#>   Agrostol   Trifrepe   Elymrepe   Anthodor   Juncarti   Trifprat   Salirepe 
#> 0.49226546 0.53434466 0.57564661 0.61543243 0.65341300 0.68921695 0.72432408 
#> 
## IGNORE_RDIFF_BEGIN
summary(sim)
#> 
#> Contrast: SF_BF 
#> 
#>          average      sd   ratio     ava     avb cumsum    p  
#> Agrostol 0.06137 0.03419 1.79490 4.66700 0.00000  0.098 0.06 .
#> Alopgeni 0.05267 0.03648 1.44390 4.33300 0.66700  0.182 0.06 .
#> Lolipere 0.04812 0.03945 1.21980 3.00000 6.00000  0.260 0.46  
#> Trifrepe 0.04630 0.02553 1.81380 1.33300 4.66700  0.334 0.10 .
#> Poatriv  0.04602 0.03380 1.36150 4.66700 3.66700  0.407 0.46  
#> Scorautu 0.04370 0.02492 1.75340 1.33300 4.33300  0.477 0.04 *
#> Bromhord 0.03368 0.02586 1.30230 0.50000 2.66700  0.531 0.02 *
#> Achimill 0.03015 0.02082 1.44820 0.16700 2.33300  0.580 0.02 *
#> Planlanc 0.02859 0.02155 1.32650 0.00000 2.00000  0.625 0.51  
#> Elymrepe 0.02807 0.02978 0.94280 2.00000 1.33300  0.670 0.57  
#> Bracruta 0.02550 0.02390 1.06690 2.00000 2.00000  0.711 0.76  
#> Poaprat  0.02513 0.02397 1.04850 2.50000 4.00000  0.751 0.84  
#> Sagiproc 0.02433 0.02215 1.09830 1.83300 0.66700  0.790 0.26  
#> Bellpere 0.01986 0.01709 1.16220 0.66700 1.66700  0.822 0.08 .
#> Eleopalu 0.01861 0.04296 0.43330 1.33300 0.00000  0.852 0.88  
#> Anthodor 0.01754 0.02580 0.67980 0.00000 1.33300  0.880 0.75  
#> Juncbufo 0.01603 0.02371 0.67620 1.16700 0.00000  0.905 0.51  
#> Vicilath 0.01467 0.01331 1.10260 0.00000 1.00000  0.929 0.02 *
#> Hyporadi 0.01029 0.01520 0.67680 0.00000 0.66700  0.945 0.56  
#> Ranuflam 0.00931 0.01360 0.68450 0.66700 0.00000  0.960 0.89  
#> Juncarti 0.00698 0.01611 0.43330 0.50000 0.00000  0.972 0.95  
#> Callcusp 0.00698 0.01611 0.43330 0.50000 0.00000  0.983 0.82  
#> Rumeacet 0.00453 0.01044 0.43330 0.33300 0.00000  0.990 0.99  
#> Cirsarve 0.00398 0.00918 0.43360 0.33300 0.00000  0.996 0.40  
#> Chenalbu 0.00233 0.00537 0.43330 0.16700 0.00000  1.000 0.38  
#> Airaprae 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA  
#> Comapalu 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA  
#> Empenigr 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA  
#> Salirepe 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA  
#> Trifprat 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Contrast: SF_HF 
#> 
#>          average      sd   ratio     ava     avb cumsum    p   
#> Agrostol 0.04738 0.03127 1.51510 4.66700 1.40000  0.084 0.29   
#> Alopgeni 0.04643 0.03290 1.41150 4.33300 1.60000  0.165 0.26   
#> Lolipere 0.04199 0.02701 1.55460 3.00000 4.00000  0.239 0.84   
#> Planlanc 0.03920 0.03321 1.18040 0.00000 3.00000  0.308 0.04 * 
#> Rumeacet 0.03899 0.02737 1.42470 0.33300 3.20000  0.377 0.01 **
#> Elymrepe 0.03188 0.02955 1.07870 2.00000 2.00000  0.433 0.38   
#> Poatriv  0.02847 0.02152 1.32270 4.66700 4.80000  0.484 1.00   
#> Bracruta 0.02526 0.02104 1.20040 2.00000 2.80000  0.528 0.93   
#> Eleopalu 0.02497 0.03888 0.64240 1.33300 0.80000  0.572 0.72   
#> Poaprat  0.02393 0.01918 1.24780 2.50000 3.40000  0.614 0.97   
#> Anthodor 0.02341 0.02143 1.09230 0.00000 1.80000  0.655 0.58   
#> Sagiproc 0.02314 0.02048 1.13010 1.83300 0.80000  0.696 0.45   
#> Trifprat 0.02308 0.02343 0.98500 0.00000 1.80000  0.737 0.01 **
#> Juncarti 0.02285 0.02568 0.88990 0.50000 1.60000  0.777 0.49   
#> Trifrepe 0.02238 0.01949 1.14860 1.33300 2.80000  0.817 0.93   
#> Juncbufo 0.02164 0.02224 0.97330 1.16700 1.20000  0.855 0.22   
#> Scorautu 0.02051 0.01642 1.24890 1.33300 2.80000  0.891 0.81   
#> Achimill 0.01518 0.01139 1.33260 0.16700 1.20000  0.918 0.72   
#> Bromhord 0.01338 0.01450 0.92220 0.50000 0.80000  0.941 0.82   
#> Ranuflam 0.01066 0.01339 0.79640 0.66700 0.40000  0.960 0.85   
#> Bellpere 0.00999 0.01257 0.79480 0.66700 0.40000  0.978 0.92   
#> Callcusp 0.00662 0.01508 0.43930 0.50000 0.00000  0.989 0.90   
#> Cirsarve 0.00381 0.00867 0.43940 0.33300 0.00000  0.996 0.55   
#> Chenalbu 0.00221 0.00503 0.43930 0.16700 0.00000  1.000 0.58   
#> Airaprae 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA   
#> Comapalu 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA   
#> Empenigr 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA   
#> Hyporadi 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA   
#> Salirepe 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA   
#> Vicilath 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Contrast: SF_NM 
#> 
#>          average      sd   ratio     ava     avb cumsum    p   
#> Poatriv  0.07828 0.04095 1.91180 4.66700 0.00000  0.101 0.01 **
#> Alopgeni 0.07122 0.04696 1.51670 4.33300 0.00000  0.194 0.01 **
#> Agrostol 0.05651 0.04418 1.27920 4.66700 2.16700  0.267 0.05 * 
#> Lolipere 0.05485 0.05991 0.91550 3.00000 0.33300  0.338 0.18   
#> Eleopalu 0.04803 0.04717 1.01820 1.33300 2.16700  0.400 0.09 . 
#> Poaprat  0.04072 0.03179 1.28100 2.50000 0.66700  0.453 0.05 * 
#> Bracruta 0.04001 0.03440 1.16310 2.00000 2.83300  0.504 0.16   
#> Elymrepe 0.03560 0.03852 0.92430 2.00000 0.00000  0.551 0.09 . 
#> Scorautu 0.03247 0.03481 0.93280 1.33300 3.16700  0.593 0.14   
#> Trifrepe 0.03043 0.03163 0.96190 1.33300 1.83300  0.632 0.70   
#> Sagiproc 0.03030 0.03048 0.99430 1.83300 0.50000  0.671 0.05 * 
#> Salirepe 0.02928 0.03201 0.91440 0.00000 1.83300  0.709 0.06 . 
#> Anthodor 0.02454 0.03669 0.66880 0.00000 1.33300  0.741 0.51   
#> Callcusp 0.02276 0.02944 0.77310 0.50000 1.16700  0.770 0.06 . 
#> Ranuflam 0.02257 0.02282 0.98890 0.66700 1.33300  0.800 0.09 . 
#> Juncarti 0.02254 0.02860 0.78830 0.50000 1.16700  0.829 0.44   
#> Hyporadi 0.02011 0.03129 0.64260 0.00000 1.16700  0.855 0.26   
#> Juncbufo 0.01986 0.02903 0.68400 1.16700 0.00000  0.881 0.28   
#> Planlanc 0.01542 0.02277 0.67720 0.00000 0.83300  0.900 0.98   
#> Airaprae 0.01488 0.02188 0.68020 0.00000 0.83300  0.920 0.09 . 
#> Bellpere 0.01232 0.01592 0.77370 0.66700 0.33300  0.936 0.79   
#> Comapalu 0.01188 0.01741 0.68260 0.00000 0.66700  0.951 0.09 . 
#> Achimill 0.00929 0.01493 0.62240 0.16700 0.33300  0.963 0.96   
#> Bromhord 0.00717 0.01633 0.43910 0.50000 0.00000  0.972 0.94   
#> Rumeacet 0.00559 0.01275 0.43840 0.33300 0.00000  0.980 0.92   
#> Empenigr 0.00523 0.01200 0.43540 0.00000 0.33300  0.986 0.30   
#> Cirsarve 0.00478 0.01089 0.43910 0.33300 0.00000  0.993 0.01 **
#> Chenalbu 0.00289 0.00660 0.43820 0.16700 0.00000  0.996 0.02 * 
#> Vicilath 0.00279 0.00642 0.43450 0.00000 0.16700  1.000 0.81   
#> Trifprat 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Contrast: BF_HF 
#> 
#>          average      sd   ratio     ava     avb cumsum    p  
#> Rumeacet 0.03867 0.02606 1.48380 0.00000 3.20000  0.082 0.05 *
#> Poatriv  0.03330 0.02579 1.29110 3.66700 4.80000  0.152 0.96  
#> Planlanc 0.03185 0.01830 1.74010 2.00000 3.00000  0.219 0.45  
#> Bromhord 0.02865 0.01799 1.59260 2.66700 0.80000  0.280 0.07 .
#> Lolipere 0.02843 0.02215 1.28340 6.00000 4.00000  0.340 0.97  
#> Elymrepe 0.02782 0.02959 0.94040 1.33300 2.00000  0.398 0.50  
#> Trifrepe 0.02584 0.01656 1.56030 4.66700 2.80000  0.453 0.81  
#> Anthodor 0.02358 0.02042 1.15470 1.33300 1.80000  0.503 0.62  
#> Achimill 0.02343 0.01474 1.58930 2.33300 1.20000  0.552 0.20  
#> Bracruta 0.02273 0.01802 1.26170 2.00000 2.80000  0.600 0.96  
#> Alopgeni 0.02161 0.02308 0.93630 0.66700 1.60000  0.646 0.91  
#> Trifprat 0.02151 0.02207 0.97470 0.00000 1.80000  0.691 0.11  
#> Juncarti 0.02008 0.02555 0.78600 0.00000 1.60000  0.734 0.69  
#> Scorautu 0.01932 0.01357 1.42410 4.33300 2.80000  0.774 0.77  
#> Bellpere 0.01829 0.01486 1.23050 1.66700 0.40000  0.813 0.17  
#> Agrostol 0.01761 0.02284 0.77080 0.00000 1.40000  0.850 1.00  
#> Juncbufo 0.01500 0.02066 0.72600 0.00000 1.20000  0.882 0.56  
#> Vicilath 0.01285 0.01140 1.12740 1.00000 0.00000  0.909 0.07 .
#> Sagiproc 0.01168 0.01297 0.90080 0.66700 0.80000  0.934 0.88  
#> Eleopalu 0.01017 0.02111 0.48170 0.00000 0.80000  0.955 0.95  
#> Hyporadi 0.00895 0.01312 0.68240 0.66700 0.00000  0.974 0.69  
#> Poaprat  0.00720 0.01010 0.71330 4.00000 3.40000  0.989 1.00  
#> Ranuflam 0.00508 0.01055 0.48170 0.00000 0.40000  1.000 0.98  
#> Airaprae 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA  
#> Chenalbu 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA  
#> Cirsarve 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA  
#> Comapalu 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA  
#> Empenigr 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA  
#> Salirepe 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA  
#> Callcusp 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Contrast: BF_NM 
#> 
#>          average      sd   ratio     ava     avb cumsum    p   
#> Lolipere 0.09068 0.02644 3.42900 6.00000 0.33300  0.124 0.01 **
#> Poatriv  0.05468 0.04465 1.22500 3.66700 0.00000  0.199 0.23   
#> Poaprat  0.05251 0.01813 2.89700 4.00000 0.66700  0.271 0.02 * 
#> Trifrepe 0.05129 0.02756 1.86100 4.66700 1.83300  0.342 0.03 * 
#> Bromhord 0.03969 0.02920 1.35900 2.66700 0.00000  0.396 0.01 **
#> Bracruta 0.03572 0.02869 1.24500 2.00000 2.83300  0.445 0.35   
#> Eleopalu 0.03376 0.03573 0.94500 0.00000 2.16700  0.491 0.40   
#> Agrostol 0.03345 0.03473 0.96300 0.00000 2.16700  0.537 0.83   
#> Achimill 0.03319 0.02338 1.42000 2.33300 0.33300  0.582 0.01 **
#> Scorautu 0.03136 0.02026 1.54800 4.33300 3.16700  0.625 0.25   
#> Anthodor 0.02806 0.03295 0.85200 1.33300 1.33300  0.664 0.29   
#> Planlanc 0.02732 0.02193 1.24600 2.00000 0.83300  0.701 0.57   
#> Salirepe 0.02677 0.02927 0.91400 0.00000 1.83300  0.738 0.12   
#> Bellpere 0.02353 0.01909 1.23200 1.66700 0.33300  0.770 0.02 * 
#> Hyporadi 0.02172 0.02450 0.88600 0.66700 1.16700  0.800 0.25   
#> Ranuflam 0.02031 0.02275 0.89300 0.00000 1.33300  0.828 0.25   
#> Elymrepe 0.01999 0.02926 0.68300 1.33300 0.00000  0.855 0.79   
#> Callcusp 0.01783 0.02681 0.66500 0.00000 1.16700  0.880 0.38   
#> Juncarti 0.01769 0.02600 0.68100 0.00000 1.16700  0.904 0.69   
#> Vicilath 0.01577 0.01447 1.09000 1.00000 0.16700  0.925 0.01 **
#> Sagiproc 0.01543 0.01857 0.83100 0.66700 0.50000  0.947 0.89   
#> Airaprae 0.01341 0.01969 0.68100 0.00000 0.83300  0.965 0.29   
#> Comapalu 0.01074 0.01571 0.68400 0.00000 0.66700  0.980 0.27   
#> Alopgeni 0.01000 0.01463 0.68300 0.66700 0.00000  0.993 0.99   
#> Empenigr 0.00479 0.01105 0.43300 0.00000 0.33300  1.000 0.46   
#> Chenalbu 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA   
#> Cirsarve 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA   
#> Juncbufo 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA   
#> Rumeacet 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA   
#> Trifprat 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Contrast: HF_NM 
#> 
#>          average      sd   ratio     ava     avb cumsum    p   
#> Poatriv  0.07155 0.01368 5.23000 4.80000 0.00000  0.099 0.01 **
#> Lolipere 0.05453 0.02962 1.84100 4.00000 0.33300  0.175 0.18   
#> Rumeacet 0.04655 0.03081 1.51100 3.20000 0.00000  0.239 0.01 **
#> Poaprat  0.04175 0.01885 2.21500 3.40000 0.66700  0.297 0.02 * 
#> Planlanc 0.04163 0.02956 1.40800 3.00000 0.83300  0.355 0.01 **
#> Bracruta 0.03534 0.02010 1.75800 2.80000 2.83300  0.404 0.34   
#> Eleopalu 0.03204 0.03231 0.99200 0.80000 2.16700  0.448 0.52   
#> Agrostol 0.03192 0.02889 1.10500 1.40000 2.16700  0.492 0.94   
#> Trifrepe 0.03037 0.02287 1.32800 2.80000 1.83300  0.534 0.63   
#> Elymrepe 0.02981 0.03868 0.77100 2.00000 0.00000  0.576 0.47   
#> Anthodor 0.02872 0.02480 1.15800 1.80000 1.33300  0.615 0.24   
#> Juncarti 0.02741 0.02854 0.96100 1.60000 1.16700  0.653 0.21   
#> Trifprat 0.02584 0.02597 0.99500 1.80000 0.00000  0.689 0.01 **
#> Salirepe 0.02534 0.02729 0.92900 0.00000 1.83300  0.724 0.16   
#> Alopgeni 0.02446 0.03240 0.75500 1.60000 0.00000  0.758 0.92   
#> Scorautu 0.02070 0.01412 1.46600 2.80000 3.16700  0.787 0.82   
#> Ranuflam 0.01928 0.01994 0.96700 0.40000 1.33300  0.814 0.36   
#> Juncbufo 0.01818 0.02465 0.73800 1.20000 0.00000  0.839 0.45   
#> Hyporadi 0.01714 0.02655 0.64600 0.00000 1.16700  0.863 0.42   
#> Callcusp 0.01683 0.02490 0.67600 0.00000 1.16700  0.886 0.40   
#> Achimill 0.01656 0.01490 1.11100 1.20000 0.33300  0.909 0.69   
#> Sagiproc 0.01528 0.01653 0.92400 0.80000 0.50000  0.930 0.93   
#> Airaprae 0.01261 0.01824 0.69100 0.00000 0.83300  0.947 0.31   
#> Bromhord 0.01209 0.01517 0.79700 0.80000 0.00000  0.964 0.79   
#> Comapalu 0.01011 0.01456 0.69400 0.00000 0.66700  0.978 0.24   
#> Bellpere 0.00880 0.01373 0.64100 0.40000 0.33300  0.990 0.91   
#> Empenigr 0.00454 0.01033 0.43900 0.00000 0.33300  0.997 0.55   
#> Vicilath 0.00240 0.00546 0.43900 0.00000 0.16700  1.000 0.89   
#> Chenalbu 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA   
#> Cirsarve 0.00000 0.00000     NaN 0.00000 0.00000  1.000   NA   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Permutation: free
#> Number of permutations: 99
## IGNORE_RDIFF_END
```

# Kendall coefficient of concordance

Function `kendall.global` computes and tests the coefficient of
concordance among several judges (variables, species) through a
permutation test.

Function `kendall.post` carries out *a posteriori* tests of the
contributions of individual judges (variables, species) to the overall
concordance of their group through permutation tests.

If several groups of judges are identified in the data table,
coefficients of concordance (`kendall.global`) or a posteriori tests
(`kendall.post`) will be computed for each group separately. Use in
ecology: to identify significant species associations.

## Usage

``` r
kendall.global(Y, group, nperm = 999, mult = "holm")
kendall.post(Y, group, nperm = 999, mult = "holm")
```

## Arguments

- Y:

  Data file (data frame or matrix) containing quantitative or
  semiquantitative data. Rows are objects and columns are judges
  (variables). In community ecology, that table is often a
  site-by-species table.

- group:

  A vector defining how judges should be divided into groups. See
  example below. If groups are not explicitly defined, all judges in the
  data file will be considered as forming a single group.

- nperm:

  Number of permutations to be performed. Default is 999.

- mult:

  Correct P-values for multiple testing using the alternatives described
  in [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) and in addition
  `"sidak"` (see Details). The Bonferroni correction is overly
  conservative; it is not recommended. It is included to allow
  comparisons with the other methods.

## Details

`Y` must contain quantitative data. They will be transformed to ranks
within each column before computation of the coefficient of concordance.

The search for species associations described in Legendre (2005)
proceeds in 3 steps:

\(1\) Correlation analysis of the species. A possible method is to
compute Ward's agglomerative clustering of a matrix of correlations
among the species. In detail: (1.1) compute a Pearson or Spearman
correlation matrix (`correl.matrix`) among the species; (1.2) turn it
into a distance matrix: `mat.D = as.dist(1-correl.matrix)`; (1.3) carry
out Ward's hierarchical clustering of that matrix using `hclust`:
`clust.ward = hclust(mat.D, "ward")`; (1.4) plot the dendrogram:
`plot(clust.ward, hang=-1)`; (1.5) cut the dendrogram in two groups,
retrieve the vector of species membership:
`group.2 = cutree(clust.ward, k=2)`. (1.6) After steps 2 and 3 below,
you may have to come back and try divisions of the species into k = \\3,
4, 5, \dots\\ groups.

\(2\) Compute global tests of significance of the 2 (or more) groups
using the function `kendall.global` and the vector defining the groups.
Groups that are not globally significant must be refined or abandoned.

\(3\) Compute a posteriori tests of the contribution of individual
species to the concordance of their group using the function
`kendall.post` and the vector defining the groups. If some species have
negative values for "Spearman.mean", this means that these species
clearly do not belong to the group, hence that group is too inclusive.
Go back to (1.5) and cut the dendrogram more finely. The left and right
groups can be cut separately, independently of the levels along the
dendrogram; write your own vector of group membership if `cutree` does
not produce the desired groups.

The corrections used for multiple testing are applied to the list of
P-values (P); they take into account the number of tests (k) carried out
simultaneously (number of groups in `kendall.global`, or number of
species in `kendall.post`). The corrections are performed using function
[`p.adjust`](https://rdrr.io/r/stats/p.adjust.html); see that function
for the description of the correction methods. In addition, there is
Šidák correction which defined as \\P\_{corr} = 1 -(1 - P)^k\\.

## Value

A table containing the following information in rows. The columns
correspond to the groups of "judges" defined in vector "group". When
function `Kendall.post` is used, there are as many tables as the number
of predefined groups.

- W :

  Kendall's coefficient of concordance, W.

- F :

  F statistic. F = W\*(m-1)/(1-W) where m is the number of judges.

- Prob.F :

  Probability associated with the F statistic, computed from the F
  distribution with nu1 = n-1-(2/m) and nu2 = nu1\*(m-1); n is the
  number of objects.

- Corrected prob.F :

  Probabilities associated with F, corrected using the method selected
  in parameter `mult`. Shown only if there are more than one group.

- Chi2 :

  Friedman's chi-square statistic (Friedman 1937) used in the
  permutation test of W.

- Prob.perm :

  Permutational probabilities, uncorrected.

- Corrected prob.perm :

  Permutational probabilities corrected using the method selected in
  parameter `mult`. Shown only if there are more than one group.

- Spearman.mean :

  Mean of the Spearman correlations between the judge under test and all
  the other judges in the same group.

- W.per.species :

  Contribution of the judge under test to the overall concordance
  statistic for that group.

## References

Friedman, M. 1937. The use of ranks to avoid the assumption of normality
implicit in the analysis of variance. Journal of the American
Statistical Association 32: 675-701.

Kendall, M. G. and B. Babington Smith. 1939. The problem of m rankings.
Annals of Mathematical Statistics 10: 275-287.

Legendre, P. 2005. Species associations: the Kendall coefficient of
concordance revisited. Journal of Agricultural, Biological, and
Environmental Statistics 10: 226-245.

Legendre, P. 2009. Coefficient of concordance. In: Encyclopedia of
Research Design. SAGE Publications (in press).

Siegel, S. and N. J. Castellan, Jr. 1988. Nonparametric statistics for
the behavioral sciences. 2nd edition. McGraw-Hill, New York.

## See also

[`cor`](https://rdrr.io/r/stats/cor.html),
[`friedman.test`](https://rdrr.io/r/stats/friedman.test.html),
[`hclust`](https://rdrr.io/r/stats/hclust.html),
[`cutree`](https://rdrr.io/r/stats/cutree.html),
[`kmeans`](https://rdrr.io/r/stats/kmeans.html),
[`cascadeKM`](https://vegandevs.github.io/vegan/reference/cascadeKM.md).

## Author

F. Guillaume Blanchet, University of Alberta, and Pierre Legendre,
Université de Montréal

## Examples

``` r
data(mite)
mite.hel <- decostand(mite, "hel")

# Reproduce the results shown in Table 2 of Legendre (2005), a single group
mite.small <- mite.hel[c(4,9,14,22,31,34,45,53,61,69),c(13:15,23)]
kendall.global(mite.small, nperm=49)
#> $Concordance_analysis
#>               Group.1
#> W          0.44160305
#> F          2.37252221
#> Prob.F     0.04403791
#> Chi2      15.89770992
#> Prob.perm  0.04000000
#> 
#> attr(,"class")
#> [1] "kendall.global"
kendall.post(mite.small, mult="holm", nperm=49)
#> $A_posteriori_tests
#>                     TVEL      ONOV      SUCT   Trhypch1
#> Spearman.mean  0.3265678 0.3965503 0.4570402 -0.1681251
#> W.per.species  0.4949258 0.5474127 0.5927802  0.1239061
#> Prob           0.1000000 0.0400000 0.0200000  0.8400000
#> Corrected prob 0.2000000 0.1200000 0.0800000  0.8400000
#> 
#> $Correction.type
#> [1] "holm"
#> 
#> attr(,"class")
#> [1] "kendall.post"

# Reproduce the results shown in Tables 3 and 4 of Legendre (2005), 2 groups
group <-c(1,1,2,1,1,1,1,1,2,1,1,1,1,1,1,2,1,2,1,1,1,1,2,1,2,1,1,1,1,1,2,2,2,2,2)
kendall.global(mite.hel, group=group, nperm=49)
#> $Concordance_analysis
#>                          Group.1      Group.2
#> W                   3.097870e-01 2.911888e-01
#> F                   1.032305e+01 4.108130e+00
#> Prob.F              1.177138e-85 4.676566e-22
#> Corrected prob.F    2.354275e-85 4.676566e-22
#> Chi2                5.130073e+02 2.210123e+02
#> Prob.perm           2.000000e-02 2.000000e-02
#> Corrected prob.perm 4.000000e-02 4.000000e-02
#> 
#> $Correction.type
#> [1] "holm"
#> 
#> attr(,"class")
#> [1] "kendall.global"
kendall.post(mite.hel, group=group, mult="holm", nperm=49)
#> $A_posteriori_tests_Group
#> $A_posteriori_tests_Group[[1]]
#>                   Brachy      PHTH     RARD      SSTR   Protopl      MEGR
#> Spearman.mean  0.1851177 0.4258111 0.359058 0.2505486 0.1802160 0.2833298
#> W.per.species  0.2190711 0.4497357 0.385764 0.2817757 0.2143736 0.3131911
#> Prob           0.0200000 0.0200000 0.020000 0.0200000 0.0200000 0.0200000
#> Corrected prob 0.7000000 0.7000000 0.700000 0.7000000 0.7000000 0.7000000
#>                      MPRO      HMIN     HMIN2      NPRA      TVEL      ONOV
#> Spearman.mean  0.09248024 0.2444656 0.4138494 0.1263751 0.4177343 0.3301159
#> W.per.species  0.13029357 0.2759462 0.4382723 0.1627761 0.4419954 0.3580278
#> Prob           0.02000000 0.0200000 0.0200000 0.0200000 0.0200000 0.0200000
#> Corrected prob 0.70000000 0.7000000 0.7000000 0.7000000 0.7000000 0.7000000
#>                     SUCT Oribatl1      PWIL  Galumna1  Stgncrs2      HRUF
#> Spearman.mean  0.2185421 0.421216 0.2574779 0.4180699 0.3623428 0.1250230
#> W.per.species  0.2511028 0.445332 0.2884163 0.4423170 0.3889118 0.1614804
#> Prob           0.0200000 0.020000 0.0200000 0.0200000 0.0200000 0.0600000
#> Corrected prob 0.7000000 0.700000 0.7000000 0.7000000 0.7000000 0.7000000
#>                     PPEL      SLAT      FSET  Lepidzts  Eupelops  Miniglmn
#> Spearman.mean  0.2188216 0.3016159 0.4217606 0.2577037 0.1108022 0.2301430
#> W.per.species  0.2513707 0.3307153 0.4458539 0.2886327 0.1478521 0.2622203
#> Prob           0.0200000 0.0200000 0.0200000 0.0200000 0.0400000 0.0200000
#> Corrected prob 0.7000000 0.7000000 0.7000000 0.7000000 0.7000000 0.7000000
#> 
#> $A_posteriori_tests_Group[[2]]
#>                     HPAV      TVIE      LCIL  Ceratoz1  Trhypch1      NCOR
#> Spearman.mean  0.1222579 0.2712078 0.1906408 0.1375601 0.1342409 0.3342345
#> W.per.species  0.2020527 0.3374616 0.2642189 0.2159637 0.2129463 0.3947586
#> Prob           0.0400000 0.0200000 0.0400000 0.0400000 0.0600000 0.0200000
#> Corrected prob 0.7000000 0.7000000 0.7000000 0.7000000 0.7000000 0.7000000
#>                     LRUG     PLAG2  Ceratoz3  Oppiminu  Trimalc2
#> Spearman.mean  0.3446561 0.1833099 0.3188922 0.1764232 0.2498877
#> W.per.species  0.4042328 0.2575544 0.3808111 0.2512938 0.3180797
#> Prob           0.0200000 0.0200000 0.0200000 0.0200000 0.0200000
#> Corrected prob 0.7000000 0.7000000 0.7000000 0.7000000 0.7000000
#> 
#> 
#> $Correction.type
#> [1] "holm"
#> 
#> attr(,"class")
#> [1] "kendall.post"

# NOTE: 'nperm' argument usually needs to be larger than 49.
# It was set to this low value for demonstration purposes.
```

# Functional Diversity and Community Distances from Species Trees

Functional diversity is defined as the total branch length in a trait
dendrogram connecting all species, but excluding the unnecessary root
segments of the tree (Petchey and Gaston 2006). Tree distance is the
increase in total branch length when combining two sites.

## Usage

``` r
treedive(comm, tree, match.force = TRUE, verbose = TRUE)
treeheight(tree)
treedist(x, tree, relative = TRUE, match.force = TRUE, ...)
```

## Arguments

- comm, x:

  Community data frame or matrix.

- tree:

  A dendrogram which for `treedive` must be for species (columns).

- match.force:

  Force matching of column names in data (`comm`, `x`) and labels in
  `tree`. If `FALSE`, matching only happens when dimensions differ (with
  a warning or message). The order of data must match to the order in
  `tree` if matching by names is not done.

- verbose:

  Print diagnostic messages and warnings.

- relative:

  Use distances relative to the height of combined tree.

- ...:

  Other arguments passed to functions (ignored).

## Details

Function `treeheight` finds the sum of lengths of connecting segments in
a dendrogram produced by
[`hclust`](https://rdrr.io/r/stats/hclust.html), or other dendrogram
that can be coerced to a correct type using
[`as.hclust`](https://rdrr.io/r/stats/as.hclust.html). When applied to a
clustering of species traits, this is a measure of functional diversity
(Petchey and Gaston 2002, 2006), and when applied to phylogenetic trees
this is phylogenetic diversity.

Function `treedive` finds the `treeheight` for each site (row) of a
community matrix. The function uses a subset of dendrogram for those
species that occur in each site, and excludes the tree root if that is
not needed to connect the species (Petchey and Gaston 2006). The subset
of the dendrogram is found by first calculating
[`cophenetic`](https://rdrr.io/r/stats/cophenetic.html) distances from
the input dendrogram, then reconstructing the dendrogram for the subset
of the cophenetic distance matrix for species occurring in each site.
Diversity is 0 for one species, and `NA` for empty communities.

Function `treedist` finds the dissimilarities among trees. Pairwise
dissimilarity of two trees is found by combining species in a common
tree and seeing how much of the tree height is shared and how much is
unique. With `relative = FALSE` the dissimilarity is defined as \\2 (A
\cup B) - A - B\\, where \\A\\ and \\B\\ are heights of component trees
and \\A \cup B\\ is the height of the combined tree. With
`relative = TRUE` the dissimilarity is \\(2(A \cup B)-A-B)/(A \cup B)\\.
Although the latter formula is similar to Jaccard dissimilarity (see
[`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.md),
[`designdist`](https://vegandevs.github.io/vegan/reference/designdist.md)),
it is not in the range \\0 \ldots 1\\, since combined tree can add a new
root. When two zero-height trees are combined into a tree of above zero
height, the relative index attains its maximum value \\2\\. The
dissimilarity is zero from a combined zero-height tree.

The functions need a dendrogram of species traits or phylogenies as an
input. If species traits contain
[`factor`](https://rdrr.io/r/base/factor.html) or
[`ordered`](https://rdrr.io/r/base/factor.html) factor variables, it is
recommended to use Gower distances for mixed data (function
[`daisy`](https://rdrr.io/pkg/cluster/man/daisy.html) in package
cluster), and usually the recommended clustering method is UPGMA
(`method = "average"` in function
[`hclust`](https://rdrr.io/r/stats/hclust.html)) (Podani and Schmera
2006). Phylogenetic trees can be changed into dendrograms using function
`as.hclust.phylo` in the ape package.

It is possible to analyse the non-randomness of tree diversity using
[`oecosimu`](https://vegandevs.github.io/vegan/reference/oecosimu.md).
This needs specifying an adequate Null model, and the results will
change with this choice.

## Value

A vector of diversity values or a single tree height, or a dissimilarity
structure that inherits from [`dist`](https://rdrr.io/r/stats/dist.html)
and can be used similarly.

## References

Lozupone, C. and Knight, R. 2005. UniFrac: a new phylogenetic method for
comparing microbial communities. *Applied and Environmental
Microbiology* 71, 8228–8235.

Petchey, O.L. and Gaston, K.J. 2002. Functional diversity (FD), species
richness and community composition. *Ecology Letters* 5, 402–411.

Petchey, O.L. and Gaston, K.J. 2006. Functional diversity: back to
basics and looking forward. *Ecology Letters* 9, 741–758.

Podani J. and Schmera, D. 2006. On dendrogram-based methods of
functional diversity. *Oikos* 115, 179–185.

## Author

Jari Oksanen

## See also

Function `treedive` is similar to the phylogenetic diversity function
`pd` in the package picante, but excludes tree root if that is not
needed to connect species. Function `treedist` is similar to the
phylogenetic similarity `phylosor` in the package picante, but excludes
unneeded tree root and returns distances instead of similarities.

[`taxondive`](https://vegandevs.github.io/vegan/reference/taxondive.md)
is something very similar from another bubble.

## Examples

``` r
## There is no data set on species properties yet, and we demonstrate
## the methods using phylogenetic trees
data(dune)
data(dune.phylodis)
cl <- hclust(dune.phylodis)
treedive(dune, cl)
#> forced matching of 'tree' labels and 'comm' names
#>         1         2         3         4         5         6         7         8 
#>  384.0913  568.8791 1172.9455 1327.9317 1426.9067 1391.1628 1479.5062 1523.0792 
#>         9        10        11        12        13        14        15        16 
#> 1460.0423 1316.4832 1366.9960 1423.5582  895.1120 1457.2705 1505.9501 1187.5165 
#>        17        18        19        20 
#>  517.6920 1394.5162 1470.4671 1439.5571 
## Significance test using Null model communities.
## The current choice fixes numbers of species and picks species
## proportionally to their overall frequency
oecosimu(dune, treedive, "r1", tree = cl, verbose = FALSE)
#> Warning: nullmodel transformed 'comm' to binary data
#> oecosimu object
#> 
#> Call: oecosimu(comm = dune, nestfun = treedive, method = "r1", tree =
#> cl, verbose = FALSE)
#> 
#> nullmodel method ‘r1’ with 99 simulations
#> 
#> alternative hypothesis: statistic is less or greater than simulated values
#> 
#>    statistic       SES    mean    2.5%     50%  97.5% Pr(sim.)   
#> 1     384.09 -1.629836  856.57  466.08  948.79 1239.6     0.03 * 
#> 2     568.88 -1.966396 1184.71  658.43 1339.86 1581.5     0.01 **
#> 3    1172.95 -0.154735 1217.03  694.28 1327.55 1611.6     0.61   
#> 4    1327.93 -0.512486 1448.33  804.29 1510.76 1691.1     0.33   
#> 5    1426.91 -0.404536 1506.72  945.93 1554.94 1752.3     0.47   
#> 6    1391.16  0.283493 1316.32  752.05 1401.15 1631.7     0.93   
#> 7    1479.51  0.244044 1425.43  853.96 1485.34 1696.3     0.91   
#> 8    1523.08  0.553853 1388.18  781.05 1439.23 1710.8     0.67   
#> 9    1460.04  0.037345 1451.28  848.96 1510.66 1745.6     0.75   
#> 10   1316.48 -0.093329 1341.17  757.48 1413.33 1648.2     0.53   
#> 11   1367.00  0.881226 1108.39  599.30 1223.99 1496.1     0.55   
#> 12   1423.56  0.955746 1128.94  607.11 1255.92 1509.1     0.31   
#> 13    895.11 -1.016986 1202.74  607.01 1324.96 1547.7     0.47   
#> 14   1457.27  1.591034  967.08  477.63 1128.27 1373.6     0.01 **
#> 15   1505.95  1.587515 1037.86  580.66 1149.69 1428.0     0.01 **
#> 16   1187.52  0.490901 1040.31  587.81 1154.80 1457.1     0.93   
#> 17    517.69 -1.534848  991.08  502.98 1132.84 1374.3     0.11   
#> 18   1394.52  0.747758 1184.51  646.46 1296.56 1496.2     0.45   
#> 19   1470.47  1.287011 1111.17  617.71 1245.12 1474.6     0.09 . 
#> 20   1439.56  1.323275 1087.82  613.94 1202.96 1449.9     0.07 . 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## Phylogenetically ordered community table
dtree <- treedist(dune, cl)
tabasco(dune, hclust(dtree), cl)

## Use tree distances in distance-based RDA
dbrda(dtree ~ 1)
#> 
#> Call: dbrda(formula = dtree ~ 1)
#> 
#>               Inertia Rank RealDims
#> Total           2.183              
#> Unconstrained   2.183   19       10
#> 
#> Inertia is squared Treedist distance
#> 
#> Eigenvalues for unconstrained axes:
#>   MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
#> 1.1971 0.4546 0.2967 0.1346 0.1067 0.0912 0.0391 0.0190 
#> (Showing 8 of 19 unconstrained eigenvalues)
#> 
```

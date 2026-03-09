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
#> 1     384.09 -1.339833  796.31  383.08  696.65 1245.9     0.09 . 
#> 2     568.88 -2.094818 1206.18  660.38 1332.47 1645.0     0.01 **
#> 3    1172.95 -0.312102 1257.26  733.98 1367.06 1595.6     0.55   
#> 4    1327.93 -0.385356 1420.82  850.40 1479.84 1721.8     0.39   
#> 5    1426.91 -0.278300 1489.02  921.64 1561.33 1777.0     0.45   
#> 6    1391.16  0.452519 1274.64  757.90 1360.98 1618.1     0.85   
#> 7    1479.51  0.186999 1441.11  868.96 1451.35 1702.7     0.87   
#> 8    1523.08  0.584696 1395.03  858.16 1438.92 1676.0     0.59   
#> 9    1460.04  0.076947 1443.58  901.32 1486.70 1677.0     0.85   
#> 10   1316.48 -0.284604 1388.86  803.70 1459.24 1717.9     0.47   
#> 11   1367.00  0.731635 1141.89  604.50 1257.67 1525.4     0.55   
#> 12   1423.56  0.906518 1174.86  641.08 1260.84 1498.1     0.29   
#> 13    895.11 -1.196901 1236.63  640.95 1327.23 1588.9     0.43   
#> 14   1457.27  1.648635  964.17  525.64 1019.98 1422.8     0.03 * 
#> 15   1505.95  1.536364 1070.92  551.51 1205.23 1450.0     0.01 **
#> 16   1187.52  0.333570 1087.07  517.78 1218.54 1467.0     0.97   
#> 17    517.69 -1.609142  999.11  545.25 1131.14 1427.1     0.01 **
#> 18   1394.52  0.805302 1166.25  628.72 1273.59 1500.8     0.45   
#> 19   1470.47  1.067514 1140.77  606.23 1254.56 1553.4     0.27   
#> 20   1439.56  1.357867 1038.14  560.26 1139.66 1466.4     0.11   
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

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
#>    statistic      SES    mean    2.5%     50%  97.5% Pr(sim.)   
#> 1     384.09 -1.53146  833.94  382.13  749.33 1239.9     0.11   
#> 2     568.88 -2.35173 1215.27  699.17 1326.75 1543.1     0.01 **
#> 3    1172.95 -0.16957 1218.58  654.49 1339.36 1516.3     0.65   
#> 4    1327.93 -0.59597 1457.63  860.80 1512.14 1735.5     0.31   
#> 5    1426.91 -0.27379 1484.48  921.82 1527.83 1746.6     0.45   
#> 6    1391.16  0.36295 1288.81  719.28 1398.64 1668.1     0.95   
#> 7    1479.51  0.19241 1435.98  848.96 1480.75 1745.6     0.99   
#> 8    1523.08  0.72095 1336.23  757.48 1403.78 1652.3     0.49   
#> 9    1460.04  0.10969 1435.34  915.83 1487.24 1766.7     0.85   
#> 10   1316.48 -0.37688 1400.08  823.83 1443.35 1671.8     0.37   
#> 11   1367.00  0.87776 1113.32  589.96 1247.83 1469.5     0.41   
#> 12   1423.56  1.03638 1130.76  617.23 1225.07 1482.9     0.21   
#> 13    895.11 -1.35947 1249.79  697.97 1317.61 1616.8     0.35   
#> 14   1457.27  1.56470 1009.49  514.48 1133.64 1399.4     0.01 **
#> 15   1505.95  1.33986 1101.12  525.10 1241.76 1470.5     0.03 * 
#> 16   1187.52  0.39241 1072.53  576.13 1191.28 1470.2     0.99   
#> 17    517.69 -1.67315  991.83  522.92 1101.36 1431.0     0.07 . 
#> 18   1394.52  0.80734 1161.52  607.40 1279.75 1497.8     0.39   
#> 19   1470.47  1.03168 1165.67  598.28 1278.00 1494.6     0.13   
#> 20   1439.56  1.28896 1051.70  585.05 1149.04 1468.8     0.15   
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

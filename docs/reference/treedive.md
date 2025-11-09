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
#> 1     384.09 -1.20343  762.44  292.05  641.23 1234.6     0.19   
#> 2     568.88 -2.14160 1200.61  668.36 1308.30 1575.8     0.01 **
#> 3    1172.95 -0.17294 1220.10  687.90 1340.97 1522.5     0.61   
#> 4    1327.93 -0.55635 1454.27  877.61 1523.67 1731.3     0.35   
#> 5    1426.91 -0.52545 1522.71  939.30 1543.93 1747.2     0.39   
#> 6    1391.16  0.26981 1322.16  709.93 1396.37 1640.8     0.97   
#> 7    1479.51  0.23512 1422.22  817.41 1478.64 1739.9     0.99   
#> 8    1523.08  0.65935 1362.25  769.51 1431.59 1695.4     0.45   
#> 9    1460.04  0.03721 1451.58  877.82 1498.39 1718.4     0.67   
#> 10   1316.48 -0.51575 1421.66  901.14 1477.60 1721.0     0.39   
#> 11   1367.00  0.78034 1119.77  571.70 1250.34 1513.8     0.51   
#> 12   1423.56  0.90420 1147.77  616.72 1265.66 1550.7     0.31   
#> 13    895.11 -1.07510 1196.58  694.74 1333.58 1509.2     0.57   
#> 14   1457.27  1.49958 1010.87  504.93 1103.63 1399.1     0.01 **
#> 15   1505.95  1.35194 1113.07  528.78 1222.09 1472.0     0.05 * 
#> 16   1187.52  0.46553 1048.62  617.74 1164.59 1481.0     0.93   
#> 17    517.69 -1.31647  923.18  511.03  800.80 1371.1     0.09 . 
#> 18   1394.52  0.81412 1154.76  588.21 1265.29 1528.5     0.41   
#> 19   1470.47  1.14152 1145.23  619.37 1252.48 1491.5     0.13   
#> 20   1439.56  1.24146 1047.22  552.77 1166.60 1478.2     0.17   
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

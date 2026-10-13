### test spantree

data(dune)
d <- vegdist(dune)
expect_silent(tr <- spantree(d))
expect_silent(spandepth(tr))
## similar to single linkage clustering
cl <- hclust(d, "single")
expect_equal(cophenetic(as.hclust(tr)), cophenetic(cl))
## cophenetic(tr) is different: it defines links between points,
## hclust links between clusters.
expect_true(all(cophenetic(tr) >= cophenetic(cl)))
## in cluster analysis, only n-1 cophenetic distances are among points
## and match spantree distances
expect_equal(sum(cophenetic(tr) == cophenetic(cl)), tr$n - 1)

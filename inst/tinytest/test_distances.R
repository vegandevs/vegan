### Tests for distance functions. Mainly vegdist, but also betadiver
### and designdist.

### vegdist methods
dmethod <- c("manhattan", "euclidean", "canberra", "bray",
             "kulczynski", "gower", "morisita", "horn", "mountford",
             "jaccard", "raup", "binomial", "chao", "altGower", "cao",
             "mahalanobis", "clark", "chisq", "chord", "hellinger",
             "aitchison", "robust.aitchison")

### check that all functions can be run
data(dune)
dists <- list()
## 'pseudocount' is only needed and used in "aitchison"
expect_silent(for (m in dmethod)
                  dists[[m]] <- vegdist(dune, m, pseudocount=1))
## check 'maxdist' argument
expect_silent(maxdist <- sapply(dists, attr, 'maxdist'))
dmax <- sapply(dists, max)
expect_true(all(maxdist >= dmax, na.rm=T))

### similar test for betadiver (which allows defining distance by number)
bdists <- list()
expect_silent(for(i in 1:24) bdists[[i]] <- betadiver(dune, i))
expect_silent(maxdist <- sapply(bdists, attr, 'maxdist'))
dmax <- sapply(bdists, max)
## betadiver has some similarities with maxdist=0
expect_true(all(maxdist == 0 | maxdist >= dmax, na.rm=T))

## crosscheck against stats::dist
expect_equivalent(dists$manhattan, dist(dune, "manhattan"))
expect_equivalent(dists$euclidean, dist(dune, "euclidean"))
expect_equivalent(dists$canberra, dist(dune, "canberra") / ncol(dune))
expect_equivalent(vegdist(dune, "jaccard", binary = TRUE),
                  dist(dune, "binary"))
## euclidean distances of transformed data
expect_equivalent(dists$chisq, dist(decostand(dune, "chi.square")))
expect_equivalent(dists$chord, dist(decostand(dune, "norm")))
expect_equivalent(dists$hellinger, dist(decostand(dune, "hellinger")))

## vegdist methods that can be defined in designdist
expect_equivalent(dists$manhattan, designdist(dune, "A+B-2*J", "minimum"))
expect_equivalent(dists$euclidean, designdist(dune, "sqrt(A+B-2*J)",
                                              "quadratic"))
expect_equivalent(dists$bray, designdist(dune, "(A+B-2*J)/(A+B)", "minimum"))
expect_equivalent(dists$kulczynski, designdist(dune, "1-J/2/A-J/2/B",
                                               "minimum"))
expect_equivalent(dists$gower, designdist(decostand(dune, "range"),
                                          "(A+B-2*J)/P", "minimum"))
expect_equivalent(dists$jaccard, designdist(dune, "(A+B-2*J)/(A+B-J)",
                                            "minimum"))
expect_equivalent(dists$raup, designdist(dune, "1-phyper(J-1,A,P-A,B)",
                                         "binary"))
chao <- suppressWarnings(chaodist(dune, "1-U*V/(U+V-U*V)"))
chao[is.na(chao)] <- 1 # contains NaN if some U==0 & V==0
expect_equivalent(dists$chao, chao)
expect_equivalent(dists$chisq, designdist(decostand(dune, "chi.square"),
                                         "sqrt(A+B-2*J)", "quadratic"))
expect_equivalent(dists$chord, designdist(decostand(dune, "norm"),
                                         "sqrt(A+B-2*J)", "quadratic"))
expect_equivalent(dists$hellinger, designdist(decostand(dune, "hellinger"),
                                              "sqrt(A+B-2*J)", "quadratic"))
## not in vegdist, but uses src/vegdist.c
expect_equivalent(no.shared(dune), designdist(dune, "J==0"))

### designdist variants: many ways of computing Bray-Curtis
bray <- vegdist(dune, "bray", binary = TRUE)
expect_equivalent(designdist(dune, "(A+B-2*J)/(A+B)"), bray)
expect_equivalent(designdist(dune, "(b+c)/(2*a+b+c)", abcd = TRUE), bray)
expect_equivalent(designdist(dune, "gamma/alpha - 1", alphagamma = TRUE), bray)

### crosscheck betadiver/vegdist
expect_equivalent(bdists[[1]], bray)
expect_equivalent(bdists[[2]], bray)
expect_equivalent(bdists[[8]], bray)
expect_equivalent(bdists[[9]], bray)
expect_equivalent(bdists[[20]], bray)
jaccard <- vegdist(dune, "jaccard", binary = TRUE)
expect_equivalent(bdists[[15]], jaccard)
expect_equivalent(bdists[[16]], jaccard)



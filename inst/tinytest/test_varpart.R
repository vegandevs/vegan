### test varpart: cannot be really verified but we just see that they
### run with different input arguments. Euclidean/metric models are
### checked for consinstency, though. RsquareAdj is also checked
### separately (implicitly in varpart).

data(dune, dune.env)

## Chi-square distances, two components
expect_silent(varpart(dune, ~ Moisture, ~ Management, data = dune.env,
                      chisquare = TRUE, permutations = 10))

## Three components, metric (Euclidean) methods are equivalent.
expect_silent(m <- varpart(dune, ~ Moisture, ~ Management, ~ A1,
                           data = dune.env, transfo = "hellinger"))
d <- vegdist(dune, "hellinger")
expect_silent(m1 <- varpart(d, ~ Moisture, ~ Management, ~ A1,
                            data = dune.env))
## item 'ordination' varies with model and must be dropped from
## comparison, but otherwise two previous should be equal.
items <- names(m$part)[-8]
expect_equal(m$part[items], m1$part[items])
## distance-based, additive constant
expect_silent(varpart(vegdist(dune), ~ Moisture, ~ Management,
                      data = dune.env, add = TRUE))

## Four components
expect_silent(varpart(dune, ~ A1, ~ Moisture, ~ Management, ~ Use,
                      data = dune.env, scale = TRUE))

## RsquareAdj
## default vs. rda method
m <- rda(dune ~ Moisture, dune.env)
expect_equal(RsquareAdj(m$CCA$tot.chi / m$tot.chi, nobs(m), m$CCA$qrank),
             RsquareAdj(m)$adj.r.squared)
## univariate rda & lm are consistent
ml <- lm(diversity(dune) ~ Moisture, dune.env)
m <- rda(diversity(dune) ~ Moisture, dune.env)
expect_equal(RsquareAdj(m), RsquareAdj(ml))

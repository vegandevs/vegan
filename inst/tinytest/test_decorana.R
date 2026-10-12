### test decorana

data(dune)

## orthogonal CA: verify against eigenvector analysis
expect_silent(m <- decorana(dune, ira = 1))
expect_silent(m0 <- ca(dune))
expect_equivalent(c(eigenvals(m)), eigenvals(m0)[1:4])
expect_silent(sco <- scores(m, display = "sites"))
expect_silent(sco0 <- scores(m0, choices = 1:4, display = "sites",
                             scaling = "sites", hill = TRUE))
expect_silent(pro <- summary(procrustes(sco, sco0, choices = 1:4)))
## Procrustes ss=0, same axes but allowing axis reversal, residuals=0
expect_equal(pro$ss, 0)
expect_equal(abs(pro$rotation), diag(4),
             tolerance = 10*sqrt(.Machine$double.eps)) # diff 3.5e-8
expect_equivalent(pro$resid, numeric(length(pro$resid)))

## predict method
expect_equal(predict(m), predict(m0, rank = 4))

## weights
expect_equivalent(weights(m), weights(m0) * m0$grand.total)

## detrended CA: cannot be verified, but see that it runs
expect_silent(md <- decorana(dune))
expect_equivalent(eigenvals(md, "decorana")[1], eigenvals(m0)[1])
expect_silent(decorana(dune, iweigh = 1))

## beforeafter: approximate sqrt transformation, max diff will be
## sqrt(2) approximated as 4/3
x <- t(seq(1, 20)) # needs a matrix
expect_equal(max(sqrt(x) - vegan:::beforeafter(x, (1:5)^2, 1:5)),
             sqrt(2) - 4/3)
expect_equivalent(
    c(eigenvals(decorana(dune, ira=1, before=1:9, after=sqrt(1:9)))),
    eigenvals(ca(sqrt(dune)))[1:4])

## tolerance method
expect_equivalent(unclass(tolerance(m, which = "species", useN2 = FALSE)),
                  unclass(tolerance(m0, useN2 = FALSE, choices = 1:4,
                                    scaling="sites", hill = TRUE)),
                  tolerance = 10*sqrt(.Machine$double.eps)) # was 1.85e-8

## bstick
expect_equivalent(bstick(m), bstick(m0)[1:4])

### decostand tests
### NB. "clr", "rclr", "alr" are tested in tests/decostand-tests.R

data(dune)

## total
expect_silent(x <- decostand(dune, "total"))
expect_equivalent(rowSums(x), rep(1, nrow(x)))
expect_equal(decobackstand(x), dune)

## max
expect_silent(x <- decostand(dune, "max"))
expect_equivalent(apply(x, 2, max), rep(1, ncol(x)))
expect_equal(decobackstand(x), dune)

## frequency
expect_silent(x <- decostand(dune, "frequency"))
expect_equivalent(colSums(x), colSums(dune > 0))
expect_equal(decobackstand(x), dune)

## normalize
expect_silent(x <- decostand(dune, "normalize"))
expect_equivalent(rowSums(x^2), rep(1, nrow(dune)))
expect_equal(decobackstand(x), dune)

## range
expect_silent(x <- decostand(dune, "range"))
expect_equivalent(apply(x, 2, max), rep(1, ncol(dune)))
expect_equal(decobackstand(x), dune)

## rank
expect_silent(x <- decostand(dune, "rank"))
expect_true(all(rowSums(dune > 0) >= apply(x, 1, max)))
expect_error(decobackstand(x))

## rrank
expect_silent(x <- decostand(dune, "rrank"))
expect_true(all(apply(x, 1, max) <= rep(1, nrow(dune))))
expect_error(decobackstand(x))

## standardize
expect_silent(x <- decostand(dune, "standardize"))
expect_equivalent(colMeans(x), rep(0, ncol(dune)))
expect_equivalent(apply(x, 2, sd), rep(1, ncol(dune)))
expect_equal(decobackstand(x), dune)

## pa
expect_silent(x <- decostand(dune, "pa"))
expect_true(all(dune[x==1] > 0))
expect_true(all(dune[x==0] == 0))
expect_error(decobackstand(x))

## chi.square
expect_silent(x <- decostand(dune, "chi.square"))
ca <- ca(dune)
m <- wcmdscale(vegdist(x, "eucl"), w = rowSums(dune)/sum(dune), eig=TRUE)
expect_equivalent(eigenvals(m), eigenvals(ca))
expect_equal(decobackstand(x), dune)

## hellinger
expect_silent(x <- decostand(dune, "hellinger"))
expect_equivalent(rowSums(x^2), rep(1, nrow(dune)))
expect_equal(decobackstand(x), dune)

## log
expect_silent(x <- decostand(dune, "log"))
expect_equal(x[dune > 0], log2(dune[dune > 0]) + 1)
expect_equal(decobackstand(x), dune)

## wisconsin
expect_silent(x <- wisconsin(dune))
expect_equivalent(rowSums(x), rep(1, nrow(x)))
expect_equal(decobackstand(x), dune)

### decostand tests
### NB. "clr", "rclr", "alr" are tested in tests/decostand-tests.R

data(dune)

## helpers to check that standardization maintains rank order (by
## MARGIN) of original data and maintains original absences (zeros)

keepranks <- function(current, target, MARGIN)
    all(apply(current, MARGIN, rank) == apply(target, MARGIN, rank))

keepzeros <- function(current, target)
    all((current == 0) == (target == 0))

## total
expect_silent(x <- decostand(dune, "total"))
expect_equivalent(rowSums(x), rep(1, nrow(x)))
expect_true(keepranks(x, dune, 1))
expect_true(keepzeros(x, dune))
expect_equal(decobackstand(x), dune)

## max
expect_silent(x <- decostand(dune, "max"))
expect_equivalent(apply(x, 2, max), rep(1, ncol(x)))
expect_true(keepranks(x, dune, 2))
expect_true(keepzeros(x, dune))
expect_equal(decobackstand(x), dune)

## frequency
expect_silent(x <- decostand(dune, "frequency"))
expect_equivalent(colSums(x), colSums(dune > 0))
expect_true(keepranks(x, dune, 2))
expect_true(keepzeros(x, dune))
expect_equal(decobackstand(x), dune)

## normalize
expect_silent(x <- decostand(dune, "normalize"))
expect_equivalent(rowSums(x^2), rep(1, nrow(dune)))
expect_true(keepranks(x, dune, 1))
expect_true(keepzeros(x, dune))
expect_equal(decobackstand(x), dune)

## range
expect_silent(x <- decostand(dune, "range"))
expect_equivalent(apply(x, 2, max), rep(1, ncol(dune)))
expect_true(keepranks(x, dune, 2))
expect_true(keepzeros(x, dune))
expect_equal(decobackstand(x), dune)

## rank: original data cannot be reconstructed
expect_silent(x <- decostand(dune, "rank"))
expect_true(all(rowSums(dune > 0) >= apply(x, 1, max)))
expect_true(keepranks(x, dune, 1))
expect_true(keepzeros(x, dune))
expect_error(decobackstand(x)) # ranks do not retain interval data

## rrank: original data cannot be reconstructed
expect_silent(x <- decostand(dune, "rrank"))
expect_true(all(apply(x, 1, max) <= rep(1, nrow(dune))))
expect_true(keepranks(x, dune, 1))
expect_true(keepzeros(x, dune))
expect_error(decobackstand(x)) # ranks do not retain interval data

## standardize
expect_silent(x <- decostand(dune, "standardize"))
expect_equivalent(colMeans(x), rep(0, ncol(dune))) # new zero
expect_equivalent(apply(x, 2, sd), rep(1, ncol(dune)))
expect_true(keepranks(x, dune, 2))
expect_false(keepzeros(x, dune)) # does not keep zeros: mean is new zero
expect_equal(decobackstand(x), dune)

## pa: original data cannot be reconstructed
expect_silent(x <- decostand(dune, "pa"))
expect_true(all(dune[x==1] > 0))
expect_true(all(dune[x==0] == 0)) # same as keepzeros()
expect_false(keepranks(x, dune, 1)) # does not keep quantitative data
expect_error(decobackstand(x)) # quantitative data cannot be found from 0/1

## chi.square: does not (except rarely) maintain rank order of original
## data due to standardization by both the row and column totals.
expect_silent(x <- decostand(dune, "chi.square"))
expect_true(keepzeros(x, dune)) # keeps original zeros
## weighted cmdscale == correspondence analysis
ca <- ca(dune)
m <- wcmdscale(vegdist(x, "eucl"), w = rowSums(dune)/sum(dune), eig=TRUE)
expect_equivalent(eigenvals(m), eigenvals(ca))
expect_equal(decobackstand(x), dune)

## hellinger
expect_silent(x <- decostand(dune, "hellinger"))
expect_equivalent(rowSums(x^2), rep(1, nrow(dune)))
expect_true(keepranks(x, dune, 1))
expect_true(keepzeros(x, dune))
expect_equal(decobackstand(x), dune)

## log
expect_silent(x <- decostand(dune, "log"))
expect_equal(x[dune > 0], log2(dune[dune > 0]) + 1)
expect_true(all(rank(x) == rank(dune))) # either way as transforms each element
expect_true(keepzeros(x, dune))
expect_equal(decobackstand(x), dune)

## Following clr, rclr (and alr) are more completely tested in
## tests/decostand-tests.R and tests/aitchison-tests.R. Here only
## similar simple tests as for other decostand methods.

## clr: needs pseudocount, will replace input zeros with negative
## constant for each row
expect_silent(x <- decostand(dune, "clr", pseudocount = 0.1))
expect_true(keepranks(x, dune, 1)) ## keeps ranks by row
expect_true(all((x < 0) == (dune == 0))) ## replace 0 with negative value
expect_equal(decobackstand(x), dune)

## rclr: default with impute=TRUE will use matrix completion to
## replace input zeros, and consequently it won't retain original
## zeros. Neither ranks, and input data cannot be reconstructed. With
## impute=FALSE, input zeros will be NA, ranks are not kept, but input
## data can be reconstructed.
expect_silent(x <- decostand(dune, "rclr", impute = FALSE))
expect_equal(decobackstand(x), dune)

### alr: does not keep ranks, zeros and cannot be reconstructed - not
### tested here either.

## wisconsin, not in decostand, but a wrapper of two decostand
## transformations. Does not (except rarely) maintain rank order of
## original data due to standardization by both columns and rows.
expect_silent(x <- wisconsin(dune))
expect_equivalent(rowSums(x), rep(1, nrow(x)))
expect_true(keepzeros(x, dune))
expect_equal(decobackstand(x), dune)

### test procrustes

## Identical configurations have zero residual sum of squares. Rounding
## in the sum-of-squares calculation could make it negative and cause
## summary() to warn when taking its square root.
x <- 0.3 * matrix(c(0, 1, 2, 3, 0, 2, 1, 4), ncol = 2)
for (symmetric in c(FALSE, TRUE)) {
    for (scale in c(FALSE, TRUE)) {
        fit <- procrustes(x, x, scale = scale, symmetric = symmetric)
        expect_true(fit$ss >= 0)
        expect_equal(fit$ss, 0)
        expect_silent(summary(fit))
    }
}

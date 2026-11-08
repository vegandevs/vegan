### test envfit, vectorfit and factorfit against direct regression and
### class means of ordination scores

data(varespec, varechem, dune, dune.env)

## vectors in unweighted ordination: r2 and arrow direction from lm
ord <- rda(varespec)
X <- scores(ord, display = "sites", choices = 1:2)
expect_silent(fit <- envfit(ord, varechem, permutations = 0))
expect_null(fit$vectors$pvals)
for (v in names(varechem)) {
    m <- lm(varechem[[v]] ~ X)
    expect_equal(unname(fit$vectors$r[v]), summary(m)$r.squared)
    b <- coef(m)[-1]
    expect_equal(unname(fit$vectors$arrows[v, ]), unname(b / sqrt(sum(b^2))))
}
expect_equal(rowSums(fit$vectors$arrows^2), rep(1, ncol(varechem)),
             check.attributes = FALSE)

## weighted ordination: arrow directions from lm with cca site weights
ord <- cca(varespec)
X <- scores(ord, display = "sites", choices = 1:2)
w <- weights(ord)
fit <- envfit(ord, varechem, permutations = 0)
for (v in c("N", "Al", "pH")) {
    m <- lm(varechem[[v]] ~ X, weights = w)
    b <- coef(m)[-1]
    expect_equal(unname(fit$vectors$arrows[v, ]), unname(b / sqrt(sum(b^2))))
    ## r2 is cor after w-centring as in vectorfit, which re-centres without
    ## weights and so differs slightly from summary(m)$r.squared
    ## (0.21957 vs 0.21979 for N). This pins current behaviour.
    Xw <- sqrt(w) * sweep(as.matrix(X), 2, colSums(X * w) / sum(w))
    yw <- sqrt(w) * (varechem[[v]] - weighted.mean(varechem[[v]], w))
    m0 <- lm(yw ~ Xw - 1)
    expect_equal(unname(fit$vectors$r[v]), cor(fitted(m0), yw)^2)
}

## factors: centroids are weighted class means, r2 is the between-class
## share of weighted variation
ord <- cca(dune)
X <- scores(ord, display = "sites", choices = 1:2)
w <- weights(ord)
fit <- envfit(ord ~ Moisture + A1, dune.env, permutations = 0)
cnt <- sapply(split(seq_len(nrow(X)), dune.env$Moisture),
              function(i) colSums(X[i, ] * w[i]) / sum(w[i]))
expect_equal(unname(fit$factors$centroids), unname(t(cnt)))
Xc <- sweep(X, 2, colSums(X * w) / sum(w))
Xb <- t(cnt)[as.integer(dune.env$Moisture), ] - matrix(colSums(X * w) / sum(w),
                                                       nrow(X), 2, byrow = TRUE)
expect_equal(unname(fit$factors$r), sum(Xb^2 * w) / sum(Xc^2 * w))
expect_equal(labels(fit),
             list(vectors = "A1",
                  factors = paste0("Moisture", c(1, 2, 4, 5))))

## formula and default interfaces agree
fit2 <- envfit(ord, dune.env[, c("Moisture", "A1")], permutations = 0)
expect_equal(fit2$vectors, fit$vectors)
expect_equal(fit2$factors, fit$factors)

## permutation P-values: same permutations give the same result,
## and they are in (0, 1]
set.seed(4711)
perm <- shuffleSet(nrow(dune), 99)
f1 <- envfit(ord ~ Moisture + A1, dune.env, permutations = perm)
f2 <- envfit(ord ~ Moisture + A1, dune.env, permutations = perm)
expect_equal(f1$vectors$pvals, f2$vectors$pvals)
expect_equal(f1$factors$pvals, f2$factors$pvals)
p <- c(f1$vectors$pvals, f1$factors$pvals)
expect_true(all(p > 0 & p <= 1))
expect_equal(f1$vectors$r, fit$vectors$r)
EPS <- sqrt(.Machine$double.eps) # as in vectorfit and factorfit
## P-value counts the observed statistic among permutations
rp <- apply(perm, 1, function(i)
    envfit(ord, dune.env[i, "A1", drop = FALSE], permutations = 0)$vectors$r)
expect_equal(unname(f1$vectors$pvals),
             (sum(rp >= fit$vectors$r - EPS) + 1) / (nrow(perm) + 1))
rp <- apply(perm, 1, function(i)
    envfit(ord, dune.env[i, "Moisture", drop = FALSE], permutations = 0)$factors$r)
expect_equal(unname(f1$factors$pvals),
             (sum(rp >= fit$factors$r - EPS) + 1) / (nrow(perm) + 1))

## missing values
env <- dune.env
env$A1[3] <- NA
expect_error(envfit(ord ~ A1, env, permutations = 0), "missing values")
fit <- envfit(ord ~ A1 + Moisture, env, permutations = 0, na.rm = TRUE)
expect_equal(as.vector(fit$na.action), 3L)
## same as dropping the row and its weight
X <- scores(ord, display = "sites", choices = 1:2)
fit2 <- envfit(X[-3, ], env[-3, c("A1", "Moisture")], permutations = 0,
               w = weights(ord)[-3])
expect_equal(fit$vectors$r, fit2$vectors$r)
expect_equal(fit$vectors$arrows, fit2$vectors$arrows)
expect_equal(fit$factors$r, fit2$factors$r)
expect_equal(fit$factors$centroids, fit2$factors$centroids)

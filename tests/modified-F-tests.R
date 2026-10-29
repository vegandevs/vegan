## Modified pseudo-F allowing heterogeneous dispersions
suppressPackageStartupMessages(library("vegan"))

tol <- 1e-8
expect_error <- function(expr, pattern) {
    msg <- tryCatch({
        force(expr)
        NULL
    }, error = conditionMessage)
    stopifnot(!is.null(msg), grepl(pattern, msg))
}

## Unequal group sizes and a fixed permutation matrix: compare the complete
## calculation with a direct hand calculation in Euclidean geometry.
group <- factor(rep(letters[1:3], c(3, 4, 5)))
Y <- cbind(
    c(1, 2, 4, 3, 6, 8, 9, 2, 5, 7, 10, 13),
    c(3, 1, 5, 4, 8, 5, 10, 7, 6, 11, 9, 15),
    c(2, 5, 3, 7, 4, 9, 6, 8, 12, 10, 14, 11)
)
perms <- rbind(
    rev(seq_len(nrow(Y))),
    c(2:nrow(Y), 1),
    c(4:12, 1:3),
    c(1, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 12)
)

manual_F2 <- function(response, group) {
    ## rda() expresses inertia as variance, so its working response has
    ## the usual sqrt(N - 1) divisor.
    response <- scale(response, scale = FALSE) / sqrt(nrow(response) - 1)
    X <- scale(model.matrix(~ group)[, -1, drop = FALSE], scale = FALSE)
    H <- X %*% solve(crossprod(X), t(X))
    fit <- H %*% response
    residual <- response - fit
    ni <- as.vector(table(group))
    Vi <- vapply(split(seq_len(nrow(response)), group),
                 function(i) sum(residual[i, , drop = FALSE]^2), numeric(1)) /
        (ni - 1)
    component <- (1 - ni / length(group)) * Vi
    denominator <- sum(component)
    den.df <- denominator^2 / sum(component^2 / (ni - 1))
    c(F = sum(fit^2) / denominator, denominator = denominator,
      den.df = den.df, Vi)
}

rmod <- rda(Y ~ group)
ptest <- permutest(rmod, permutations = perms, statistic = "modified-F",
                  groups = group, parallel = 1)
ptest.direct <- permutest(rmod, permutations = perms, model = "direct",
                         statistic = "modified-F", groups = group,
                         parallel = 1)
hand <- manual_F2(Y, group)
stopifnot(
    identical(ptest$statistic, "modified-F"),
    all.equal(ptest$F.perm, ptest.direct$F.perm, tolerance = tol),
    all.equal(unname(ptest$F.0), unname(hand["F"]), tolerance = tol),
    all.equal(ptest$den.0, unname(hand["denominator"]), tolerance = tol),
    all.equal(ptest$den.df, unname(hand["den.df"]), tolerance = tol),
    all.equal(unname(ptest$group.dispersion), unname(hand[levels(group)]),
              tolerance = tol),
    all.equal(drop(ptest$F.perm),
              apply(perms, 1, function(i) manual_F2(Y[i, ], group)["F"]),
              tolerance = tol),
    identical(dim(permustats(ptest)$permutations), dim(ptest$F.perm))
)

## Balanced designs reduce exactly to the ordinary pseudo-F, both observed
## and under every permutation.
balanced <- factor(rep(LETTERS[1:3], each = 4))
bmod <- rda(Y ~ balanced)
b1 <- permutest(bmod, permutations = perms, groups = balanced,
                statistic = "pseudo-F", parallel = 1)
b2 <- permutest(bmod, permutations = perms, groups = balanced,
                statistic = "modified-F", parallel = 1)
stopifnot(
    !"statistic" %in% names(b1),
    all.equal(unname(b1$F.0), unname(b2$F.0), tolerance = tol),
    all.equal(drop(b1$F.perm), drop(b2$F.perm), tolerance = tol)
)

## Euclidean adonis2, dbrda and unscaled rda agree.
dmod <- dbrda(dist(Y) ~ group)
dtest <- permutest(dmod, permutations = perms, statistic = "modified-F",
                  groups = group, parallel = 1)
atest <- adonis2(Y ~ group, method = "euclidean", permutations = perms,
                 statistic = "modified-F", parallel = 1)
stopifnot(
    identical(attr(atest, "statistic"), "modified-F"),
    all(c("denominator", "denominator.df", "groups", "group.size",
          "group.dispersion") %in% names(attributes(atest))),
    all.equal(unname(ptest$F.0), unname(dtest$F.0), tolerance = tol),
    all.equal(drop(ptest$F.perm), drop(dtest$F.perm), tolerance = tol),
    all.equal(unname(ptest$F.0), unname(atest$F[1]), tolerance = tol),
    all.equal(drop(ptest$F.perm), drop(attr(atest, "F.perm")),
              tolerance = tol)
)

## Scaled rda agrees with distance methods on the correspondingly scaled
## response.
Ys <- scale(Y)
rsmod <- rda(Y ~ group, scale = TRUE)
dstest <- permutest(dbrda(dist(Ys) ~ group), permutations = perms,
                   statistic = "modified-F", groups = group, parallel = 1)
astest <- adonis2(Ys ~ group, method = "euclidean", permutations = perms,
                  statistic = "modified-F", parallel = 1)
rstest <- permutest(rsmod, permutations = perms, statistic = "modified-F",
                   groups = group, parallel = 1)
stopifnot(
    all.equal(unname(rstest$F.0), unname(dstest$F.0), tolerance = tol),
    all.equal(drop(rstest$F.perm), drop(dstest$F.perm), tolerance = tol),
    all.equal(unname(rstest$F.0), unname(astest$F[1]), tolerance = tol)
)

## A saturated factorial omnibus model equals its combined one-way cell
## model. Deliberately use unequal cell sizes.
A <- factor(rep(c("a1", "a2", "a1", "a2"), c(2, 3, 3, 4)))
B <- factor(rep(c("b1", "b1", "b2", "b2"), c(2, 3, 3, 4)))
cell <- interaction(A, B, drop = TRUE)
fac <- permutest(rda(Y ~ A * B), permutations = perms,
                 statistic = "modified-F", groups = cell, parallel = 1)
one <- permutest(rda(Y ~ cell), permutations = perms,
                 statistic = "modified-F", groups = cell, parallel = 1)
stopifnot(
    all.equal(unname(fac$F.0), unname(one$F.0), tolerance = tol),
    all.equal(drop(fac$F.perm), drop(one$F.perm), tolerance = tol)
)

## Non-Euclidean dissimilarities, restricted permutations, and serial versus
## forked parallel execution.
data(dune, dune.env)
blocks <- factor(rep(seq_len(5), each = 4))
set.seed(194)
restricted <- permute::shuffleSet(nrow(dune), 11,
                                  control = permute::how(blocks = blocks))
bray1 <- adonis2(dune ~ Management, dune.env, permutations = restricted,
                 statistic = "modified-F", parallel = 1)
bray2 <- adonis2(dune ~ Management, dune.env, permutations = restricted,
                 statistic = "modified-F", parallel = 2)
bray.dist <- as.matrix(vegdist(dune))
bray.Vi <- vapply(split(seq_len(nrow(dune)), dune.env$Management),
                  function(i) sum(bray.dist[i, i]^2) /
                      (2 * length(i) * (length(i) - 1)), numeric(1))
stopifnot(
    is.finite(bray1$F[1]),
    all.equal(unname(attr(bray1, "group.dispersion")), unname(bray.Vi),
              tolerance = tol),
    all.equal(bray1$F[1], bray2$F[1], tolerance = tol),
    all.equal(attr(bray1, "F.perm"), attr(bray2, "F.perm"), tolerance = tol)
)

## Explicit and inferred groups receive the model's NA omission.
gna <- group
gna[2] <- NA
na.data <- data.frame(gna = gna)
na.perms <- matrix(c(2:(nrow(Y) - 1L), 1L), nrow = 1L)
na.inferred <- adonis2(Y ~ gna, na.data, permutations = na.perms,
                       na.action = na.omit, method = "euclidean",
                       statistic = "modified-F")
na.explicit <- adonis2(Y ~ gna, na.data, permutations = na.perms,
                       na.action = na.omit, method = "euclidean",
                       statistic = "modified-F", groups = gna)
na.direct <- permutest(rda(Y[-2, ] ~ droplevels(gna[-2])),
                      permutations = na.perms, statistic = "modified-F",
                      groups = droplevels(gna[-2]))
stopifnot(
    all.equal(na.inferred$F[1], unname(na.direct$F.0), tolerance = tol),
    all.equal(na.inferred$F[1], na.explicit$F[1], tolerance = tol)
)

## Explicit groups can be named only in the supplied data frame.
data.only <- data.frame(data.group = group)
data.group.test <- adonis2(Y ~ data.group, data.only, permutations = perms,
                           method = "euclidean", statistic = "modified-F",
                           groups = data.group)
stopifnot(all.equal(data.group.test$F[1], unname(ptest$F.0), tolerance = tol))

## Published benchmark from Anderson et al. (2017), Section 4. For Areas 2
## and 3, presence/absence data and Jaccard dissimilarities gave ordinary
## pseudo-F = 6.53, modified F2 = 6.84, and P = 0.0001 (9999 permutations).
data("norwegianShelf", package = "vegan")
stopifnot(
    identical(dim(norwegianShelf$comm), c(101L, 809L)),
    identical(unname(as.integer(table(norwegianShelf$area))),
              c(16L, 21L, 25L, 19L, 20L)),
    !anyNA(norwegianShelf$comm),
    all(norwegianShelf$comm >= 0L)
)
take23 <- norwegianShelf$area %in% c("2", "3")
norway23 <- 1L * (norwegianShelf$comm[take23, , drop = FALSE] > 0L)
area23 <- droplevels(norwegianShelf$area[take23])
set.seed(491)
paper.modified <- adonis2(norway23 ~ area23, method = "jaccard",
                          binary = TRUE, permutations = 9999,
                          statistic = "modified-F", parallel = 1)
paper.ordinary <- adonis2(norway23 ~ area23, method = "jaccard",
                          binary = TRUE,
                          permutations = matrix(seq_len(nrow(norway23)),
                                                nrow = 1L),
                          parallel = 1)

## Check the group dispersions directly from equation (2), independently of
## the residualized-Gower implementation.
paper.dist <- as.matrix(vegdist(norway23, method = "jaccard", binary = TRUE))
paper.Vi <- vapply(split(seq_len(nrow(norway23)), area23), function(i) {
    sum(paper.dist[i, i]^2) / (2 * length(i) * (length(i) - 1))
}, numeric(1))
stopifnot(
    all.equal(unname(paper.modified$F[1]), 6.838793122991717,
              tolerance = tol),
    identical(round(paper.modified$F[1], 2), 6.84),
    identical(paper.modified[1, "Pr(>F)"], 0.0001),
    all.equal(unname(paper.ordinary$F[1]), 6.527746539283111,
              tolerance = tol),
    identical(round(paper.ordinary$F[1], 2), 6.53),
    all.equal(unname(attr(paper.modified, "group.dispersion")),
              unname(paper.Vi), tolerance = tol),
    all.equal(unname(attr(paper.modified, "denominator")),
              0.257539641107903, tolerance = tol),
    all.equal(unname(attr(paper.modified, "denominator.df")),
              43.6590318640301, tolerance = tol)
)

## Validation of the published method's scope.
singleton <- factor(c("single", as.character(group[-1])))
wrong.groups <- factor(c(rep("x", 4), rep("y", 4), rep("z", 4)))
z <- seq_len(nrow(Y))
two.groups <- factor(rep(c("lo", "hi"), each = nrow(Y) / 2))
two.numeric <- as.numeric(two.groups)
expect_error(
    permutest(rda(Y ~ singleton), permutations = perms,
              statistic = "modified-F", groups = singleton),
    "at least two observations"
)
expect_error(
    permutest(rmod, permutations = perms, statistic = "modified-F",
              groups = wrong.groups),
    "span the same space"
)
expect_error(
    permutest(rda(Y ~ A + B), permutations = perms,
              statistic = "modified-F", groups = cell),
    "omnibus one-way"
)
expect_error(
    permutest(rda(Y ~ two.numeric), permutations = perms,
              statistic = "modified-F", groups = two.groups),
    "constraint terms must be factors"
)
expect_error(
    permutest(rda(Y ~ group + Condition(z)), permutations = perms,
              statistic = "modified-F", groups = group),
    "not available with Conditions"
)
expect_error(
    anova(rmod, permutations = perms, by = "terms",
          statistic = "modified-F", groups = group),
    "'by' must be NULL"
)
expect_error(
    anova(rmod, test = "F", statistic = "modified-F", groups = group),
    "requires a permutation test"
)
expect_error(
    permutest(rmod, permutations = perms, first = TRUE,
              statistic = "modified-F", groups = group),
    "'first' must be FALSE"
)
expect_error(
    permutest(rmod, permutations = perms, model = "full",
              statistic = "modified-F", groups = group),
    "model = .*full"
)
expect_error(
    anova(rmod, rmod, permutations = perms, statistic = "modified-F",
          groups = group),
    "cannot compare multiple models"
)
expect_error(
    adonis2(Y ~ A * B, permutations = perms, method = "euclidean",
            statistic = "modified-F"),
    "must be supplied for a factorial"
)
expect_error(
    permutest(cca(Y ~ group), permutations = perms,
              statistic = "modified-F", groups = group),
    "not available for cca"
)
expect_error(
    permutest(capscale(dist(Y) ~ group), permutations = perms,
              statistic = "modified-F", groups = group),
    "not available for capscale"
)

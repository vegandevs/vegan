### permutation comparisons for constrained ordination

data(dune, dune.env)

set.seed(4711)
per <- shuffleSet(nrow(dune), 49)

### A one-response RDA comparison reduces to the ordinary partial F test
set.seed(1)
udat <- data.frame(group = factor(rep(LETTERS[1:3], each = 10)),
                   x = rnorm(30))
uy <- with(udat, 2 * (group == "B") + x + rnorm(30))
urda <- rda(uy ~ group + x, data = udat)
ucmp <- ordicomparisons(urda, "group", comparison = "B - A = 0",
                    permutations = 9)
ulm <- lm(uy ~ group + x, data = udat)
expect_equal(ucmp$F, coef(summary(ulm))["groupB", "t value"]^2,
             tolerance = 1e-8,
             info = "one-response RDA comparison matches partial F test")

### Pairwise RDA matches an explicit nested parameterization
mod <- rda(dune ~ Management, data = dune.env)
cmp <- ordicomparisons(mod, "Management", comparison = "HF - BF = 0",
                   permutations = per)

denv <- dune.env
denv$Management_BF_HF <- factor(ifelse(
    denv$Management %in% c("BF", "HF"), "BF_HF",
    as.character(denv$Management)))
denv$BF_vs_HF <- with(denv, ifelse(Management == "BF", -1,
                            ifelse(Management == "HF", 1, 0)))
restricted <- rda(dune ~ Management_BF_HF, data = denv)
full <- rda(dune ~ Management_BF_HF + BF_vs_HF, data = denv)
nested <- anova(restricted, full, permutations = per)
expect_equal(cmp$Inertia, nested[[4L]][2], tolerance = 1e-8,
             info = "comparison inertia matches nested RDA")
expect_equal(cmp$F, nested[[5L]][2], tolerance = 1e-8,
             info = "comparison pseudo-F matches nested RDA")

### Comparisons are invariant to fitted contrast coding
old <- options(contrasts = c("contr.sum", "contr.poly"))
sum.mod <- rda(dune ~ Management, data = dune.env)
options(old)
sum.cmp <- ordicomparisons(sum.mod, "Management", comparison = "HF - BF = 0",
                       permutations = per)
expect_equal(cmp$Inertia, sum.cmp$Inertia, tolerance = 1e-8,
             info = "comparison inertia is coding invariant")
expect_equal(cmp$F, sum.cmp$F, tolerance = 1e-8,
             info = "comparison pseudo-F is coding invariant")

### Profile rows expose the population being compared
profile.mod <- rda(dune ~ Management + A1, data = dune.env)
sc <- ordicontrasts(profile.mod, "Management")
named.sc <- ordicontrasts(profile.mod, variable = "Management")
expect_inherits(sc, "vegan_contrasts")
expect_identical(sc$Profile, named.sc$Profile,
                 info = "named and positional variable arguments agree")
expect_identical(sc$Profile, paste0("b", 1:4))
expect_identical(as.character(sc$Management), levels(dune.env$Management))
expect_true(all(sc$N.reference == nrow(dune)))
expect_false(any(sc$Extrapolated))
expect_identical(attr(sc, "averaged"), "A1")
expect_match(capture.output(print(sc))[4L], "Reference population:")

### Profile identifiers agree with factor-level names
bcmp <- ordicomparisons(sc, comparison = "b2 - b1 = 0", permutations = per)
named.bcmp <- ordicomparisons(profile.mod, "Management",
                          comparison = "HF - BF = 0", permutations = per)
named.arg.bcmp <- ordicomparisons(profile.mod, variable = "Management",
                              comparison = "HF - BF = 0",
                              permutations = per)
expect_equal(named.bcmp$Inertia, bcmp$Inertia, tolerance = 1e-8)
expect_equal(named.bcmp$F, bcmp$F, tolerance = 1e-8)
expect_equal(named.bcmp$F, named.arg.bcmp$F, tolerance = 1e-8)
expect_identical(attr(bcmp, "contrasts"), sc,
                 info = "comparisons retain the profiles defining the test")
expect_identical(rownames(attr(sc, "profile.matrix")), sc$Profile,
                 info = "profile identifiers label the retained design rows")
expect_identical(attr(bcmp, "comparison.matrices")$profiles,
                 attr(sc, "profile.matrix"))
expect_equal(drop(attr(bcmp, "comparison.matrices")$comparisons[[1L]]),
             drop(attr(sc, "profile.matrix")[2L, ] -
                  attr(sc, "profile.matrix")[1L, ]),
             info = "displayed profile rows define the tested difference")

### Explicit and mean-or-mode reference populations
within.sc <- ordicontrasts(profile.mod, "Management",
                                  newdata = data.frame(A1 = 5))
outside.sc <- ordicontrasts(profile.mod, "Management",
                                   newdata = data.frame(A1 = 12.5))
mean.sc <- ordicontrasts(profile.mod, "Management", newdata = "mean")
mean.explicit <- ordicontrasts(profile.mod, "Management",
                                      newdata = data.frame(
                                          A1 = mean(dune.env$A1)))
expect_false(any(within.sc$Extrapolated))
expect_true(all(outside.sc$Extrapolated))
expect_identical(unname(attr(within.sc, "extrapolation")), FALSE)
expect_identical(unname(attr(outside.sc, "extrapolation")), TRUE)
expect_true(all(mean.sc$N.reference == 1L))
expect_identical(attr(mean.sc, "averaged"), character())
expect_match(attr(mean.sc, "population"), "mean-or-mode reference profile")
expect_equal(attr(mean.sc, "profile.matrix"),
             attr(mean.explicit, "profile.matrix"), tolerance = 1e-8,
             info = "mean profile matches explicit numeric mean")
expect_error(ordicontrasts(mod, "Use"),
             info = "non-model focal variable is rejected")
expect_error(ordicomparisons(sc[1:2, ], permutations = per),
             info = "edited profile table is rejected")

### Built-in generators, equations, joint tests, and p-value adjustment
pw <- ordicomparisons(mod, "Management", permutations = per, adjust = "holm")
expect_equal(nrow(pw), choose(nlevels(dune.env$Management), 2))
expect_true("Pr(>F).adj" %in% names(pw))
expect_true(all(pw$N.Perm == nrow(per)))
expect_inherits(attr(pw, "model"), "rda")

ctrl <- ordicomparisons(mod, "Management", comparison = "trt.vs.ctrl",
                    permutations = per)
ctrl.explicit <- ordicomparisons(mod, "Management", comparison = "trt.vs.ctrl",
                             reference = "BF", permutations = per)
expect_equal(nrow(ctrl), nlevels(dune.env$Management) - 1L)
expect_identical(ctrl$Comparison, ctrl.explicit$Comparison)
expect_equal(ctrl$F, ctrl.explicit$F, tolerance = 1e-8,
             info = "default control is the first fitted factor level")

for (recipe in c("sequential", "consecutive", "reference", "deviation",
                 "effect", "helmert", "poly")) {
    ans <- ordicomparisons(mod, "Management", comparison = recipe,
                       permutations = per)
    expect_true(nrow(ans) >= 3L, info = paste("recipe", recipe, "runs"))
}

custom <- ordicomparisons(mod, "Management",
                      comparison = c("HF - BF = 0",
                                     "(NM + SF) / 2 - BF = 0"),
                      permutations = per)
expect_equal(nrow(custom), 2L)
named <- ordicomparisons(mod, "Management",
                     comparison = c(one = "HF - BF = 0",
                                    two = "NM - BF = 0"),
                     permutations = per)
expect_equal(nrow(named), 2L)
profile.custom <- ordicomparisons(ordicontrasts(mod, "Management"),
                               comparison = c("b2 - b1 = 0",
                                   "(b3 + b4) / 2 - b1 = 0"),
                               permutations = per)
expect_equal(nrow(profile.custom), 2L)
expect_equal(custom$F, profile.custom$F, tolerance = 1e-8,
             info = "profile and factor-level expressions agree")

general <- ordicomparisons(sc, comparison = "b1 + b2 - b3 - b4 = 0",
                       permutations = per)
expect_equal(nrow(general), 1L)
expect_equal(general$Df, 1L)
joint <- ordicomparisons(mod, "Management", comparison = "trt.vs.ctrl",
                     reference = "BF", joint = TRUE, permutations = per)
expect_equal(nrow(joint), 1L)
expect_equal(joint$Df, 3L)

### Interactions, explicit moderator values, and grouped comparisons
set.seed(10)
idat <- expand.grid(trt = factor(c("A", "B", "C")),
                    season = factor(c("spring", "autumn")), rep = 1:4)
Y <- matrix(rpois(nrow(idat) * 5, 5), ncol = 5)
imod <- rda(Y ~ trt * season, data = idat)
iper <- shuffleSet(nrow(idat), 49)
iavg <- ordicomparisons(imod, "trt", permutations = iper)
isc <- ordicontrasts(imod, "trt", by = "season")
imean.by <- ordicontrasts(imod, "trt", by = "season",
                                 newdata = "mean")
iat <- ordicomparisons(imod, "trt", newdata = data.frame(
    season = factor(c("spring", "autumn"), levels = levels(idat$season))),
    permutations = iper)
iby <- ordicomparisons(imod, "trt", by = "season", permutations = iper)
iby.sc <- ordicomparisons(isc, comparison = "pairwise", permutations = iper)
idiff <- ordicomparisons(isc, comparison = "(b5 - b4) - (b2 - b1) = 0",
                     permutations = iper)
expect_equal(nrow(iavg), 3L)
expect_equal(nrow(iat), 3L)
expect_equal(nrow(iby), 6L)
expect_true("season" %in% names(iby))
expect_equal(nrow(idiff), 1L)
expect_equal(iby$F, iby.sc$F, tolerance = 1e-8,
             info = "grouped model and profile interfaces agree")
expect_equal(nrow(isc), 6L)
expect_true(all(isc$N.reference == 12L))
expect_identical(attr(isc, "averaged"), character())
expect_equal(nrow(imean.by), 6L)
expect_true(all(imean.by$N.reference == 1L))
expect_equal(nrow(attr(imean.by, "reference.data")), 2L)

ibal <- ordicomparisons(imod, "trt", newdata = "balanced", permutations = iper)
iwgt <- ordicomparisons(imod, "trt", by = "season", wts = seq_len(nrow(idat)),
                    permutations = iper)
fixed <- ordicontrasts(imod, "trt", newdata = data.frame(
    season = factor("spring", levels = levels(idat$season))))
expect_equal(nrow(ibal), 3L)
expect_equal(nrow(iwgt), 6L)
expect_identical(attr(fixed, "fixed"), "season = spring")

### A mean-or-mode profile uses the modal non-focal factor level
midx <- c(seq_len(nrow(idat)), which(idat$season == "spring")[1:3])
mdat <- idat[midx, , drop = FALSE]
mY <- Y[midx, , drop = FALSE]
mmod <- rda(mY ~ trt * season, data = mdat)
modal.sc <- ordicontrasts(mmod, "trt", newdata = "mean")
modal.explicit <- ordicontrasts(mmod, "trt", newdata = data.frame(
    season = factor("spring", levels = levels(mdat$season))))
expect_identical(as.character(attr(modal.sc, "reference.data")$season),
                 "spring")
expect_equal(attr(modal.sc, "profile.matrix"),
             attr(modal.explicit, "profile.matrix"), tolerance = 1e-8,
             info = "modal factor profile matches explicit profile")
expect_error(ordicontrasts(mmod, "trt", newdata = "median"))

### Reference data reuse the transformations estimated when fitting
set.seed(2026)
tdata <- data.frame(trt = factor(rep(c("A", "B", "C"), each = 20)),
                    x = runif(60, -2, 3),
                    season = factor(rep(c("spring", "autumn"), 30),
                                    levels = c("spring", "autumn")))
tY <- matrix(rnorm(nrow(tdata) * 4), ncol = 4)
tY[, 1] <- tY[, 1] + (tdata$trt == "B") * tdata$x
tY[, 2] <- tY[, 2] + (tdata$trt == "C") * tdata$x^2
tper <- shuffleSet(nrow(tdata), 9)
tgrid <- data.frame(x = c(-1, 0.25, 2),
                    season = factor(c("spring", "autumn", "spring"),
                                    levels = levels(tdata$season)))
tweights <- c(1, 2, 1)

for (transformation in c("poly", "scale", "square")) {
    if (transformation == "poly") {
        tbasis <- poly(tdata$x, 2)
        tformula <- tY ~ trt * poly(x, 2) + season
    } else if (transformation == "scale") {
        tbasis <- scale(tdata$x)
        tformula <- tY ~ trt * scale(x) + season
    } else {
        tbasis <- matrix(tdata$x^2, ncol = 1)
        tformula <- tY ~ trt * I(x^2) + season
    }
    basis.names <- paste0("basis", seq_len(ncol(tbasis)))
    for (j in seq_along(basis.names))
        tdata[[basis.names[j]]] <- tbasis[, j]
    fixed.formula <- reformulate(
        paste0("trt * (", paste(basis.names, collapse = " + "), ") + season"),
        response = "tY")
    transformed.fit <- rda(tformula, data = tdata)
    fixed.fit <- rda(fixed.formula, data = tdata)

    for (population in c("observed", "mean", "balanced", "custom")) {
        if (population == "observed") {
            requested <- fixed.grid <- NULL
            weights <- NULL
        } else {
            if (population == "mean") {
                requested <- "mean"
                values <- data.frame(x = mean(tdata$x),
                                     season = factor("spring",
                                         levels = levels(tdata$season)))
            } else if (population == "balanced") {
                requested <- "balanced"
                values <- data.frame(x = mean(tdata$x),
                                     season = factor(levels(tdata$season),
                                         levels = levels(tdata$season)))
            } else {
                requested <- values <- tgrid
            }
            new.basis <- if (transformation == "poly")
                predict(tbasis, values$x)
            else if (transformation == "scale")
                (matrix(values$x, ncol = 1) - attr(tbasis, "scaled:center")) /
                    attr(tbasis, "scaled:scale")
            else matrix(values$x^2, ncol = 1)
            fixed.grid <- values["season"]
            for (j in seq_along(basis.names))
                fixed.grid[[basis.names[j]]] <- new.basis[, j]
            weights <- if (population == "custom") tweights else NULL
        }
        actual <- ordicomparisons(transformed.fit, "trt",
                     comparison = "B - A = 0", newdata = requested,
                     wts = weights, permutations = tper)
        expected <- ordicomparisons(fixed.fit, "trt",
                       comparison = "B - A = 0", newdata = fixed.grid,
                       wts = weights, permutations = tper)
        expect_equal(actual$Inertia, expected$Inertia, tolerance = 1e-8,
                     info = paste(transformation, population,
                                  "inertia agrees with fixed fitted basis"))
        expect_equal(actual$F, expected$F, tolerance = 1e-8,
                     info = paste(transformation, population,
                                  "F agrees with fixed fitted basis"))
        expect_equal(actual[["Pr(>F)"]], expected[["Pr(>F)"]],
                     info = paste(transformation, population,
                                  "permutations agree with fixed fitted basis"))
    }

    ### Transformations are fitted before subset and missing-data removal
    tmissing <- tdata
    tmissing$season[c(3, 29)] <- NA
    keep <- seq_len(nrow(tmissing)) %% 5 != 0
    selected.fit <- rda(tformula, data = tmissing, subset = keep,
                         na.action = na.omit)
    selected.fixed <- rda(fixed.formula, data = tmissing, subset = keep,
                           na.action = na.omit)
    use <- keep & !is.na(tmissing$season)
    selected.mean <- mean(tmissing$x[use])
    new.basis <- if (transformation == "poly") predict(tbasis, selected.mean)
                 else if (transformation == "scale")
                     matrix((selected.mean - attr(tbasis, "scaled:center")) /
                              attr(tbasis, "scaled:scale"), ncol = 1)
                 else matrix(selected.mean^2, ncol = 1)
    selected.mode <- names(which.max(table(tmissing$season[use])))
    fixed.grid <- data.frame(season = factor(selected.mode,
                                            levels = levels(tdata$season)))
    for (j in seq_along(basis.names))
        fixed.grid[[basis.names[j]]] <- new.basis[, j]
    selected.per <- shuffleSet(sum(use), 9)
    actual <- ordicomparisons(selected.fit, "trt", newdata = "mean",
                 comparison = "B - A = 0", permutations = selected.per)
    expected <- ordicomparisons(selected.fixed, "trt", newdata = fixed.grid,
                   comparison = "B - A = 0", permutations = selected.per)
    expect_equal(actual$Inertia, expected$Inertia, tolerance = 1e-8,
                 info = paste(transformation, "preserves the fitted basis after subset and NA handling"))
    expect_equal(actual$F, expected$F, tolerance = 1e-8)
    expect_equal(actual[["Pr(>F)"]], expected[["Pr(>F)"]])
}

### Wrappers and qualified calls keep the fitted transformation parameters
wrapped.formulas <- list(
    scale = tY ~ trt * I(scale(x)) + season,
    poly = tY ~ trt * I(poly(x, 2)[, 1]) + season,
    qualified = tY ~ trt * base::scale(x) + season)
for (transformation in names(wrapped.formulas)) {
    if (transformation == "poly") {
        tbasis <- poly(tdata$x, 2)
        tdata$fixedbasis <- tbasis[, 1]
        at.mean <- predict(tbasis, mean(tdata$x))[, 1]
    } else {
        tbasis <- scale(tdata$x)
        tdata$fixedbasis <- as.vector(tbasis)
        at.mean <- (mean(tdata$x) - attr(tbasis, "scaled:center")) /
            attr(tbasis, "scaled:scale")
    }
    wrapped.fit <- rda(wrapped.formulas[[transformation]], data = tdata)
    fixed.fit <- rda(tY ~ trt * fixedbasis + season, data = tdata)
    actual <- ordicomparisons(wrapped.fit, "trt", newdata = "mean",
                 comparison = "B - A = 0", permutations = tper)
    expected <- ordicomparisons(fixed.fit, "trt", newdata = data.frame(
                    fixedbasis = as.vector(at.mean),
                    season = factor("spring", levels = levels(tdata$season))),
                    comparison = "B - A = 0", permutations = tper)
    expect_equal(actual$Inertia, expected$Inertia, tolerance = 1e-8,
                 info = paste(transformation, "wrapper retains the fitted basis"))
    expect_equal(actual$F, expected$F, tolerance = 1e-8)
    expect_equal(actual[["Pr(>F)"]], expected[["Pr(>F)"]])
}

### A data column named scale does not hide the transformation function
named.data <- tdata
named.data$scale <- seq_len(nrow(named.data))
named.fit <- rda(tY ~ trt * I(scale(x)) + season, data = named.data)
plain.fit <- rda(tY ~ trt * I(scale(x)) + season, data = tdata)
named.profiles <- ordicontrasts(named.fit, "trt", newdata = "mean")
plain.profiles <- ordicontrasts(plain.fit, "trt", newdata = "mean")
expect_equal(attr(named.profiles, "profile.matrix"),
             attr(plain.profiles, "profile.matrix"), tolerance = 1e-8,
             info = "data columns do not hide fitted transformation functions")

### A polynomial of a scaled variable retains both fitted transformations
scaled.x <- scale(tdata$x)
nested.basis <- poly(as.numeric(scaled.x), 2)
tdata$nested1 <- nested.basis[, 1]
tdata$nested2 <- nested.basis[, 2]
nested.fit <- rda(tY ~ trt * poly(scale(x), 2) + season, data = tdata)
nested.fixed <- rda(tY ~ trt * (nested1 + nested2) + season, data = tdata)
for (population in c("mean", "custom")) {
    values <- if (population == "mean")
        data.frame(x = mean(tdata$x),
                   season = factor("spring", levels = levels(tdata$season)))
    else tgrid
    scaled.values <- (values$x - attr(scaled.x, "scaled:center")) /
        attr(scaled.x, "scaled:scale")
    new.basis <- predict(nested.basis, scaled.values)
    fixed.grid <- data.frame(season = values$season,
                             nested1 = new.basis[, 1],
                             nested2 = new.basis[, 2])
    actual <- ordicomparisons(nested.fit, "trt",
                 newdata = if (population == "mean") "mean" else values,
                 comparison = "B - A = 0", permutations = tper)
    expected <- ordicomparisons(nested.fixed, "trt", newdata = fixed.grid,
                   comparison = "B - A = 0", permutations = tper)
    expect_equal(actual$Inertia, expected$Inertia, tolerance = 1e-8,
                 info = paste("polynomial of scaled x", population,
                              "retains both fitted transformations"))
    expect_equal(actual$F, expected$F, tolerance = 1e-8)
    expect_equal(actual[["Pr(>F)"]], expected[["Pr(>F)"]])
}

### Older fits can be used unless fitted transformations need recovering
old.simple <- profile.mod
old.simple$terminfo$constraint.terms <- NULL
old.simple.cmp <- ordicomparisons(old.simple, "Management",
                     comparison = "HF - BF = 0", permutations = per)
expect_equal(old.simple.cmp$Inertia, named.bcmp$Inertia, tolerance = 1e-8)
expect_equal(old.simple.cmp$F, named.bcmp$F, tolerance = 1e-8,
             info = "older fits without transformations remain usable")
old.transformed <- wrapped.fit
old.transformed$terminfo$constraint.terms <- NULL
expect_error(ordicomparisons(old.transformed, "trt", newdata = "mean",
                            permutations = tper),
             pattern = "refit the model",
             info = "older transformed fits must be refitted to retain their basis")

### Conditions and CCA use weighted chi-square geometry
ccamod <- cca(dune ~ Management + Condition(Moisture), data = dune.env)
ccacmp <- ordicomparisons(ccamod, "Management", permutations = per)
expect_equal(nrow(ccacmp), 6L)
expect_true(all(is.finite(ccacmp$F)))
cca.restricted <- cca(dune ~ Management_BF_HF + Condition(Moisture), data = denv)
cca.full <- cca(dune ~ Management_BF_HF + BF_vs_HF + Condition(Moisture),
                data = denv)
cca.nested <- anova(cca.restricted, cca.full, permutations = per)
expect_equal(ccacmp$Inertia[1], cca.nested[[4L]][2], tolerance = 1e-8,
             info = "comparison inertia agrees with explicit nested partial CCA")
expect_equal(ccacmp$F[1], cca.nested[[5L]][2], tolerance = 1e-8,
             info = "comparison F agrees with explicit nested partial CCA")
old <- options(contrasts = c("contr.sum", "contr.poly"))
cca.sum <- cca(dune ~ Management + Condition(Moisture), data = dune.env)
options(old)
cca.sum.cmp <- ordicomparisons(cca.sum, "Management", permutations = per)
expect_equal(ccacmp$Inertia, cca.sum.cmp$Inertia, tolerance = 1e-8,
             info = "partial CCA comparison inertia is coding invariant")
expect_equal(ccacmp$F, cca.sum.cmp$F, tolerance = 1e-8,
             info = "partial CCA comparison F is coding invariant")

### Euclidean dbRDA and retained adonis2 fits agree with RDA
dbmod <- dbrda(dune ~ Management, data = dune.env, distance = "euclidean")
dbcmp <- ordicomparisons(dbmod, "Management", comparison = "HF - BF = 0",
                     permutations = per)
afit <- adonis2(dune ~ Management, data = dune.env, method = "euclidean",
                return = "model")
acmp <- ordicomparisons(afit, "Management", comparison = "HF - BF = 0",
                    permutations = per)
expect_inherits(afit, "adonis2_fit")
expect_equal(cmp$F, dbcmp$F, tolerance = 1e-8,
             info = "Euclidean dbRDA agrees with RDA")
expect_equal(cmp$F, acmp$F, tolerance = 1e-8,
             info = "retained adonis2 fit agrees with RDA")

### ANOVA bridge
atab <- anova(mod, by = "contrast", variable = "Management",
              comparison = "reference", reference = "BF",
              permutations = per)
expect_inherits(atab, "anova.cca")
expect_equal(nrow(atab), 3L)
expect_error(anova(mod, by = "contrast", variable = "Management",
                   test = "parametric"),
             info = "contrast bridge rejects parametric F test")
expect_error(anova(mod, by = "contrast", variable = "Management",
                   permustat = "modified-F", groups = dune.env$Management),
             pattern = "contrast comparisons require",
             info = "modified-F is not silently ignored by contrast bridge")

### Retained fits preserve the upstream modified-F omnibus calculation
modified.table <- adonis2(dune ~ Management, data = dune.env,
                          method = "euclidean", permutations = per,
                          permustat = "modified-F")
modified.retained <- anova(afit, permutations = per,
                           permustat = "modified-F",
                           groups = dune.env$Management)
expect_equal(modified.table$F[1], modified.retained$F[1], tolerance = 1e-8)
expect_equal(modified.table[["Pr(>F)"]][1],
             modified.retained[["Pr(>F)"]][1])

### Invalid and unsafe comparisons are rejected
for (recipe in c("revpairwise", "revsequential", "revreference")) {
    expect_error(ordicomparisons(mod, "Management", comparison = recipe,
                             permutations = per),
                 info = paste("obsolete recipe", recipe, "is rejected"))
}
expect_error(ordicomparisons(mod, "Management", comparison = "HF / BF = 0",
                         permutations = per),
             info = "division by a profile is rejected")
expect_error(ordicomparisons(mod, "Management",
                         comparison = "system('false') = 0",
                         permutations = per),
             info = "function calls in expressions are rejected")
expect_error(ordicomparisons(mod, "Management", comparison = "HF - BF = 1",
                         permutations = per),
             info = "nonzero null is rejected")

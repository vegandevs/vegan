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

### Scenario rows expose the population being compared
scenario.mod <- rda(dune ~ Management + A1, data = dune.env)
sc <- ordicontrasts(scenario.mod, "Management")
named.sc <- ordicontrasts(scenario.mod, variable = "Management")
expect_inherits(sc, "vegan_scenarios")
expect_identical(sc$Scenario, named.sc$Scenario,
                 info = "named and positional variable arguments agree")
expect_identical(sc$Scenario, paste0("b", 1:4))
expect_identical(as.character(sc$Management), levels(dune.env$Management))
expect_true(all(sc$N.reference == nrow(dune)))
expect_false(any(sc$Extrapolated))
expect_identical(attr(sc, "averaged"), "A1")
expect_match(capture.output(print(sc))[4L], "Reference population:")

### Scenario identifiers agree with factor-level names
bcmp <- ordicomparisons(sc, comparison = "b2 - b1 = 0", permutations = per)
named.bcmp <- ordicomparisons(scenario.mod, "Management",
                          comparison = "HF - BF = 0", permutations = per)
named.arg.bcmp <- ordicomparisons(scenario.mod, variable = "Management",
                              comparison = "HF - BF = 0",
                              permutations = per)
expect_equal(named.bcmp$Inertia, bcmp$Inertia, tolerance = 1e-8)
expect_equal(named.bcmp$F, bcmp$F, tolerance = 1e-8)
expect_equal(named.bcmp$F, named.arg.bcmp$F, tolerance = 1e-8)

### Explicit and mean-or-mode reference populations
within.sc <- ordicontrasts(scenario.mod, "Management",
                                  newdata = data.frame(A1 = 5))
outside.sc <- ordicontrasts(scenario.mod, "Management",
                                   newdata = data.frame(A1 = 12.5))
mean.sc <- ordicontrasts(scenario.mod, "Management", newdata = "mean")
mean.explicit <- ordicontrasts(scenario.mod, "Management",
                                      newdata = data.frame(
                                          A1 = mean(dune.env$A1)))
expect_false(any(within.sc$Extrapolated))
expect_true(all(outside.sc$Extrapolated))
expect_identical(unname(attr(within.sc, "extrapolation")), FALSE)
expect_identical(unname(attr(outside.sc, "extrapolation")), TRUE)
expect_true(all(mean.sc$N.reference == 1L))
expect_identical(attr(mean.sc, "averaged"), character())
expect_match(attr(mean.sc, "population"), "mean-or-mode reference profile")
expect_equal(attr(mean.sc, "scenario.matrix"),
             attr(mean.explicit, "scenario.matrix"), tolerance = 1e-8,
             info = "mean profile matches explicit numeric mean")
expect_error(ordicontrasts(mod, "Use"),
             info = "non-model focal variable is rejected")
expect_error(ordicomparisons(sc[1:2, ], permutations = per),
             info = "edited scenario table is rejected")

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
scenario.custom <- ordicomparisons(ordicontrasts(mod, "Management"),
                               comparison = c("b2 - b1 = 0",
                                   "(b3 + b4) / 2 - b1 = 0"),
                               permutations = per)
expect_equal(nrow(scenario.custom), 2L)
expect_equal(custom$F, scenario.custom$F, tolerance = 1e-8,
             info = "scenario and factor-level expressions agree")

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
             info = "grouped model and scenario interfaces agree")
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
expect_equal(attr(modal.sc, "scenario.matrix"),
             attr(modal.explicit, "scenario.matrix"), tolerance = 1e-8,
             info = "modal factor profile matches explicit profile")
expect_error(ordicontrasts(mmod, "trt", newdata = "median"))

### Conditions and CCA use weighted chi-square geometry
ccamod <- cca(dune ~ Management + Condition(Moisture), data = dune.env)
ccacmp <- ordicomparisons(ccamod, "Management", permutations = per)
expect_equal(nrow(ccacmp), 6L)
expect_true(all(is.finite(ccacmp$F)))

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
             info = "division by a scenario is rejected")
expect_error(ordicomparisons(mod, "Management",
                         comparison = "system('false') = 0",
                         permutations = per),
             info = "function calls in expressions are rejected")
expect_error(ordicomparisons(mod, "Management", comparison = "HF - BF = 1",
                         permutations = per),
             info = "nonzero null is rejected")

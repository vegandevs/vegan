## permutation comparisons for constrained ordination

library(vegan)
library(permute)

data(dune)
data(dune.env)

set.seed(4711)
per <- shuffleSet(nrow(dune), 49)

## A one-response RDA comparison reduces to the ordinary partial F test.
set.seed(1)
udat <- data.frame(group = factor(rep(LETTERS[1:3], each = 10)),
                   x = rnorm(30))
uy <- with(udat, 2 * (group == "B") + x + rnorm(30))
urda <- rda(uy ~ group + x, data = udat)
ucmp <- comparisons(urda, "group", comparison = "B - A = 0",
                    permutations = 9)
ulm <- lm(uy ~ group + x, data = udat)
stopifnot(all.equal(ucmp$F,
                    coef(summary(ulm))["groupB", "t value"]^2,
                    tolerance = 1e-8))

## Pairwise RDA comparison: observed statistic agrees with an explicit
## nested parameterization that retains every treatment level.
mod <- rda(dune ~ Management, data = dune.env)
cmp <- comparisons(mod, "Management", comparison = "HF - BF = 0",
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
stopifnot(all.equal(cmp$Inertia, nested[[4L]][2], tolerance = 1e-8),
          all.equal(cmp$F, nested[[5L]][2], tolerance = 1e-8))

## The scientific comparison is invariant to the fitted contrast coding.
old <- options(contrasts = c("contr.sum", "contr.poly"))
sum.mod <- rda(dune ~ Management, data = dune.env)
options(old)
sum.cmp <- comparisons(sum.mod, "Management", comparison = "HF - BF = 0",
                       permutations = per)
stopifnot(all.equal(cmp$Inertia, sum.cmp$Inertia, tolerance = 1e-8),
          all.equal(cmp$F, sum.cmp$F, tolerance = 1e-8))

## Scenario rows expose the population being compared. Positional names refer
## to these rows, and are equivalent to unambiguous factor-level names.
scenario.mod <- rda(dune ~ Management + A1, data = dune.env)
sc <- model_scenarios(scenario.mod, "Management")
stopifnot(inherits(sc, "vegan_scenarios"),
          identical(sc$Scenario, paste0("b", 1:4)),
          identical(as.character(sc$Management), levels(dune.env$Management)),
          all(sc$N.reference == nrow(dune)), !any(sc$Extrapolated),
          identical(attr(sc, "averaged"), "A1"),
          grepl("Reference population:", capture.output(print(sc))[4L]))
bcmp <- comparisons(sc, comparison = "b2 - b1 = 0", permutations = per)
named.bcmp <- comparisons(scenario.mod, "Management",
                          comparison = "HF - BF = 0", permutations = per)
stopifnot(all.equal(named.bcmp$Inertia, bcmp$Inertia, tolerance = 1e-8),
          all.equal(named.bcmp$F, bcmp$F, tolerance = 1e-8))
within.sc <- model_scenarios(scenario.mod, "Management",
                             newdata = data.frame(A1 = 5))
outside.sc <- model_scenarios(scenario.mod, "Management",
                              newdata = data.frame(A1 = 12.5))
stopifnot(!any(within.sc$Extrapolated), all(outside.sc$Extrapolated),
          identical(unname(attr(within.sc, "extrapolation")), FALSE),
          identical(unname(attr(outside.sc, "extrapolation")), TRUE))
stopifnot(inherits(try(model_scenarios(mod, "Use"), silent = TRUE),
                     "try-error"),
          inherits(try(comparisons(sc[1:2, ], permutations = per),
                       silent = TRUE), "try-error"))

## Built-in generators, safe equations, scenario expressions, joint tests,
## and p.adjust.
pw <- comparisons(mod, "Management", permutations = per, adjust = "holm")
stopifnot(nrow(pw) == choose(nlevels(dune.env$Management), 2),
          "Pr(>F).adj" %in% names(pw), all(pw$N.Perm == nrow(per)),
          inherits(attr(pw, "model"), "rda"))
ctrl <- comparisons(mod, "Management", comparison = "trt.vs.ctrl",
                    permutations = per)
ctrl.explicit <- comparisons(mod, "Management", comparison = "trt.vs.ctrl",
                             reference = "BF", permutations = per)
stopifnot(nrow(ctrl) == nlevels(dune.env$Management) - 1L,
          identical(ctrl$Comparison, ctrl.explicit$Comparison),
          all.equal(ctrl$F, ctrl.explicit$F, tolerance = 1e-8))
for (recipe in c("sequential", "consecutive", "reference", "deviation",
                 "effect", "helmert", "poly"))
    stopifnot(nrow(comparisons(mod, "Management", comparison = recipe,
                              permutations = per)) >= 3L)
custom <- comparisons(mod, "Management",
                      comparison = c("HF - BF = 0",
                                     "(NM + SF) / 2 - BF = 0"),
                      permutations = per)
stopifnot(nrow(custom) == 2L)
named <- comparisons(mod, "Management",
                     comparison = c(one = "HF - BF = 0",
                                    two = "NM - BF = 0"),
                     permutations = per)
stopifnot(nrow(named) == 2L)
scenario.custom <- comparisons(model_scenarios(mod, "Management"),
                               comparison = c("b2 - b1 = 0",
                                   "(b3 + b4) / 2 - b1 = 0"),
                               permutations = per)
stopifnot(nrow(scenario.custom) == 2L,
          all.equal(custom$F, scenario.custom$F, tolerance = 1e-8))
general <- comparisons(sc, comparison = "b1 + b2 - b3 - b4 = 0",
                       permutations = per)
stopifnot(nrow(general) == 1L, general$Df == 1L)
joint <- comparisons(mod, "Management", comparison = "trt.vs.ctrl",
                     reference = "BF", joint = TRUE, permutations = per)
stopifnot(nrow(joint) == 1L, joint$Df == 3L)

## Interactions: empirical averages, explicit moderator values, and by groups.
set.seed(10)
idat <- expand.grid(trt = factor(c("A", "B", "C")),
                    season = factor(c("spring", "autumn")), rep = 1:4)
Y <- matrix(rpois(nrow(idat) * 5, 5), ncol = 5)
imod <- rda(Y ~ trt * season, data = idat)
iper <- shuffleSet(nrow(idat), 49)
iavg <- comparisons(imod, "trt", permutations = iper)
isc <- model_scenarios(imod, "trt", by = "season")
iat <- comparisons(imod, "trt", newdata = data.frame(
    season = factor(c("spring", "autumn"), levels = levels(idat$season))),
    permutations = iper)
iby <- comparisons(imod, "trt", by = "season",
                   permutations = iper)
iby.sc <- comparisons(isc, comparison = "pairwise", permutations = iper)
idiff <- comparisons(isc, comparison = "(b5 - b4) - (b2 - b1) = 0",
                     permutations = iper)
stopifnot(nrow(iavg) == 3L, nrow(iat) == 3L, nrow(iby) == 6L,
          "season" %in% names(iby), nrow(idiff) == 1L,
          all.equal(iby$F, iby.sc$F, tolerance = 1e-8),
          nrow(isc) == 6L, all(isc$N.reference == 12L),
          identical(attr(isc, "averaged"), character()))
ibal <- comparisons(imod, "trt", newdata = "balanced", permutations = iper)
iwgt <- comparisons(imod, "trt", by = "season", wts = seq_len(nrow(idat)),
                    permutations = iper)
fixed <- model_scenarios(imod, "trt", newdata = data.frame(
    season = factor("spring", levels = levels(idat$season))))
stopifnot(nrow(ibal) == 3L, nrow(iwgt) == 6L,
          identical(attr(fixed, "fixed"), "season = spring"))

## Conditions and CCA use the same compiler with weighted chi-square geometry.
ccamod <- cca(dune ~ Management + Condition(Moisture), data = dune.env)
ccacmp <- comparisons(ccamod, "Management", permutations = per)
stopifnot(nrow(ccacmp) == 6L, all(is.finite(ccacmp$F)))

## Euclidean dbRDA and retained adonis2 fits agree with RDA pseudo-F values.
dbmod <- dbrda(dune ~ Management, data = dune.env, distance = "euclidean")
dbcmp <- comparisons(dbmod, "Management", comparison = "HF - BF = 0",
                     permutations = per)
afit <- adonis2(dune ~ Management, data = dune.env, method = "euclidean",
                return = "model")
acmp <- comparisons(afit, "Management", comparison = "HF - BF = 0",
                    permutations = per)
stopifnot(inherits(afit, "adonis2_fit"),
          all.equal(cmp$F, dbcmp$F, tolerance = 1e-8),
          all.equal(cmp$F, acmp$F, tolerance = 1e-8))

## ANOVA bridge.
atab <- anova(mod, by = "contrast", variables = "Management",
              comparison = "reference", reference = "BF",
              permutations = per)
stopifnot(inherits(atab, "anova.cca"), nrow(atab) == 3L)

## Invalid and unsafe comparisons are rejected.
stopifnot(all(vapply(c("revpairwise", "revsequential", "revreference"),
                     function(recipe)
                         inherits(try(comparisons(
                             mod, "Management", comparison = recipe,
                             permutations = per), silent = TRUE),
                                  "try-error"), logical(1))),
          inherits(try(comparisons(mod, "Management",
                                  comparison = "HF / BF = 0",
                                  permutations = per), silent = TRUE),
                   "try-error"),
          inherits(try(comparisons(mod, "Management",
                                  comparison = "system('false') = 0",
                                  permutations = per), silent = TRUE),
                   "try-error"),
          inherits(try(comparisons(mod, "Management",
                                  comparison = "HF - BF = 1",
                                  permutations = per), silent = TRUE),
                   "try-error"))

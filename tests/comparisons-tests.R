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

## Built-in generators, safe equations, composition, joint tests and p.adjust.
pw <- comparisons(mod, "Management", permutations = per, adjust = "holm")
stopifnot(nrow(pw) == choose(nlevels(dune.env$Management), 2),
          "Pr(>F).adj" %in% names(pw), all(pw$N.Perm == nrow(per)),
          inherits(attr(pw, "model"), "rda"))
ref <- comparisons(mod, "Management", comparison = "trt.vs.ctrl",
                   reference = "BF", permutations = per)
stopifnot(nrow(ref) == nlevels(dune.env$Management) - 1L)
for (recipe in c("revpairwise", "sequential", "consecutive",
                 "revsequential", "reference", "revreference",
                 "deviation", "effect", "helmert", "poly"))
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
                     hypothesis = "two - one = 0", permutations = per)
stopifnot(nrow(named) == 1L)
composed <- comparisons(mod, "Management", comparison = "pairwise",
                        hypothesis = "b2 - b1 = 0", permutations = per)
stopifnot(nrow(composed) == 1L, composed$Df == 1L)
general <- comparisons(mod, "Management", comparison = "pairwise",
                       hypothesis = "b1 + b2 - b3 = 0", permutations = per)
stopifnot(nrow(general) == 1L)
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
iat <- comparisons(imod, "trt", newdata = data.frame(
    season = factor(c("spring", "autumn"), levels = levels(idat$season))),
    permutations = iper)
iby <- comparisons(imod, "trt", by = "season",
                   permutations = iper)
idiff <- comparisons(imod, "trt", by = "season",
                     hypothesis = "b4 - b1 = 0", permutations = iper)
stopifnot(nrow(iavg) == 3L, nrow(iat) == 3L, nrow(iby) == 6L,
          "season" %in% names(iby), nrow(idiff) == 1L)
ibal <- comparisons(imod, "trt", newdata = "balanced", permutations = iper)
iwgt <- comparisons(imod, "trt", by = "season", wts = seq_len(nrow(idat)),
                    permutations = iper)
stopifnot(nrow(ibal) == 3L, nrow(iwgt) == 6L)

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

## Invalid and unsafe hypotheses are rejected.
stopifnot(inherits(try(comparisons(mod, "Management",
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

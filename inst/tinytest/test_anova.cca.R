### verification tests for anova.cca cases (incl. adonis2)

data(dune, dune.env)
data(mite, mite.env) # more variable marginal totals to see cca re-weighting

## Crosscheck anova/permutest implementation. getFcore was the R
## prototype used to write C code (getF.c) for CCA re-weighting, and
## current compiled code was written after this prototype. XgetFcore
## is an alternative implementation permuting constraints. It should
## give the same result as getFcore. I considered switching getF.c to
## use this algorithm, but have not done so yet: it could be faster
## for cca, but slower for unweighted methods (rda, dbrda, adonis2).
source("R/getFcore.R")
perm <- shuffleSet(nrow(mite), 19)
mod <- cca(mite ~ Condition(Topo + Shrub) + WatrCont + SubsDens, mite.env)
expect_silent(ano <- anova(mod, permutations = perm))
proto <- getFcore(mod, perm) # returns F-values
expect_equal(drop(permustats(ano)$permutations), proto)
xproto <- XgetFcore(mod, perm)
expect_equal(drop(permustats(ano)$permutations), xproto)

perm <- shuffleSet(nrow(dune), 49)
### All anova methods run smoothly
mod0 <- rda(dune ~ Management + Moisture, dune.env)
expect_silent(anova(mod0, permutations = 49))
expect_silent(anova(mod0, by="term", permutations = 49))
expect_silent(anova(mod0, by = "margin", permutations = 49))
expect_silent(anova(mod0, by ="onedf", permutations = 49))
expect_silent(anova(mod0, by = "axis", permutations = 49))
expect_silent(anova(mod0, by = "margin", test = "parametric"))

### Euclidean methods
an0 <- anova(mod0, permutations = perm)
pstat0 <- permustats(an0)$permutations
d <- dist(dune)
ano <- anova(dbrda(d ~ Management + Moisture, dune.env), permutations=perm)
expect_equal(permustats(ano)$permutations, pstat0,
             info="rda and Euclidean dbrda match")
ano <- anova(capscale(d ~ Management + Moisture, dune.env), permutations=perm)
expect_equal(permustats(ano)$permutations, pstat0,
             info="rda and Euclidean capscale match")
ado <- adonis2(d ~ Management + Moisture, dune.env, permutations = perm)
expect_equal(permustats(ado)$permutations, pstat0,
             info="rda and Euclidean adonis2 match")

### Distance-based method with adjustments
d <- vegdist(dune)
an0 <- anova(dbrda(d ~ Management + Moisture, dune.env), permutations = perm)
ado <- adonis2(d ~ Management + Moisture, dune.env, permutations = perm)
expect_equal(permustats(an0)$permutations, permustats(ado)$permutations,
             info="dbrda and adonis2 match")
## sqrt.dist=TRUE
an0 <- anova(dbrda(d ~ Management + Moisture, dune.env, sqrt.dist=TRUE),
             permutations = perm)
ado <- adonis2(d ~ Management + Moisture, dune.env, sqrt.dist = TRUE,
               permutations = perm)
ano <- anova(capscale(d ~ Management + Moisture, dune.env, sqrt.dist=TRUE),
             permutations=perm)
expect_equal(permustats(an0)$permutations, permustats(ado)$permutations,
             info="dbrda and adonis2 match w/sqrt.dist")
expect_equal(permustats(an0)$permutations, permustats(ano)$permutations,
             info="dbrda and capscale match w/sqrt.dist")
## add = TRUE (Lingoes)
an0 <- anova(dbrda(d ~ Management + Moisture, dune.env, add=TRUE),
             permutations = perm)
ado <- adonis2(d ~ Management + Moisture, dune.env, add = TRUE,
               permutations = perm)
ano <- anova(capscale(d ~ Management + Moisture, dune.env, add=TRUE),
             permutations=perm)
expect_equal(permustats(an0)$permutations, permustats(ado)$permutations,
             info="dbrda and adonis2 match w/ add=TRUE")
expect_equal(permustats(an0)$permutations, permustats(ano)$permutations,
             info="dbrda and capscale match w/ add=TRUE")
## add = TRUE (Lingoes)
an0 <- anova(dbrda(d ~ Management + Moisture, dune.env, add=TRUE),
             permutations = perm)
ado <- adonis2(d ~ Management + Moisture, dune.env, add = TRUE,
               permutations = perm)
ano <- anova(capscale(d ~ Management + Moisture, dune.env, add=TRUE),
             permutations = perm)
expect_equal(permustats(an0)$permutations, permustats(ado)$permutations,
             info="dbrda and adonis2 match w/ add=TRUE")
expect_equal(permustats(an0)$permutations, permustats(ano)$permutations,
             info="dbrda and capscale match w/ add=TRUE")
# add = cailliez
an0 <- anova(dbrda(d ~ Management + Moisture, dune.env, add="cailliez"),
             permutations = perm)
ado <- adonis2(d ~ Management + Moisture, dune.env, add = "cailliez",
               permutations = perm)
ano <- anova(capscale(d ~ Management + Moisture, dune.env, add="cailliez"),
             permutations = perm)
expect_equal(permustats(an0)$permutations, permustats(ado)$permutations,
             info="dbrda and adonis2 match w/ add='cailliez'")
expect_equal(permustats(an0)$permutations, permustats(ano)$permutations,
             info="dbrda and capscale match w/ add='cailliez'")

### Wide, rank-deficient QR decompositions (#814)
dune6 <- dune[1:6, ]
d6 <- dist(dune6)
group <- rep(0:1, each = 3)
trend <- rep(-1:1, 2)
wide <- outer(trend, seq_len(8L))
perm6 <- rbind(c(2L, 4L, 1L, 6L, 3L, 5L),
               c(4L, 2L, 6L, 1L, 5L, 3L),
               c(6L, 5L, 3L, 4L, 2L, 1L),
               c(3L, 1L, 5L, 2L, 6L, 4L))

## Use one matrix-valued term so dropping it removes all its columns.
for (by in list(NULL, "terms", "margin", "onedf")) {
    expect_message(ado <- adonis2(d6 ~ group + wide, permutations = perm6, by = by),
                   "Some constraints or conditions were aliased")
    ref <- adonis2(d6 ~ group + trend, permutations = perm6, by = by)
    expect_equal(unname(as.matrix(ado)), unname(as.matrix(ref)))
    expect_equal(unname(permustats(ado)$permutations),
                 unname(permustats(ref)$permutations))
}

## Both the conditions and the combined model have wide QR matrices.
expect_message(ado <- adonis2(d6 ~ group + Condition(wide), permutations = perm6),
               "Some constraints or conditions were aliased")
ref <- adonis2(d6 ~ group + Condition(trend), permutations = perm6)
expect_equal(as.matrix(ado), as.matrix(ref))
expect_equal(permustats(ado)$permutations, permustats(ref)$permutations)

mod <- suppressMessages(cca(dune6 ~ group + Condition(wide)))
ano <- anova(mod, permutations = perm6)
ref <- anova(cca(dune6 ~ group + Condition(trend)), permutations = perm6)
expect_equal(as.matrix(ano), as.matrix(ref))
expect_equal(permustats(ano)$permutations, permustats(ref)$permutations)

## Undefined overall tests leave the permutation argument unevaluated.
observation <- factor(seq_len(nrow(dune6)))
expect_warning(ado <- adonis2(d6 ~ observation,
                             permutations = stop("unexpected permutation request")),
               "No residual degrees of freedom")
expect_equal(ado$Df, c(5, 0, 5))
expect_true(all(is.na(ado$F) & is.na(ado$`Pr(>F)`)))
expect_null(attr(ado, "F.perm"))

expect_warning(ado <- adonis2(d6 ~ trend + Condition(trend),
                             permutations = stop("unexpected permutation request")),
               "No estimable constraints")
expect_equal(ado$Df, c(0, 4, 5))
expect_true(all(is.na(ado$`Pr(>F)`)))
expect_null(attr(ado, "F.perm"))

expect_warning(ado <- adonis2(d6 ~ observation + Condition(observation),
                             permutations = perm6),
               "No residual degrees of freedom")
expect_equal(ado$Df, c(0, 0, 5))

## Duplicated terms have zero marginal df and the other test is unchanged.
duplicate <- trend
for (p in list(perm6, perm6[1, , drop = FALSE])) {
    expect_warning(ado <- adonis2(d6 ~ trend + duplicate + group,
                                 permutations = p, by = "margin"),
                   "Marginal tests have zero degrees of freedom for .*trend.*duplicate")
    expect_equal(ado$Df[1:2], c(0, 0))
    expect_true(all(is.na(ado$F[1:2]) & is.na(ado$`Pr(>F)`[1:2])))
    ref <- adonis2(d6 ~ trend + group, permutations = p, by = "margin")
    expect_equal(as.matrix(ado["group", ]), as.matrix(ref["group", ]))
    expect_equal(permustats(ado)$statistic, permustats(ref)$statistic["group"])
    expect_equal(permustats(ado)$permutations,
                 permustats(ref)$permutations[, 2, drop = FALSE])
}
expect_warning(ado <- adonis2(d6 ~ trend + duplicate,
                             permutations = perm6, by = "margin"),
               "Marginal tests have zero degrees of freedom")
expect_true(all(is.na(ado$F) & is.na(ado$`Pr(>F)`)))
expect_null(attr(ado, "F.perm"))

### 'by' special cases
m0 <- cca(dune ~ 1, dune.env)
m1 <- cca(dune ~ Management, dune.env)
m2 <- update(m1, . ~ . + Moisture, dune.env)
anoterm <- anova(m2, by = "term", permutations = perm)
anolist <- anova(m0, m1, m2, permutations = perm)
expect_equal(permustats(anoterm)$permutations,
             permustats(anolist)$permutations,
             info = "anova(..., by='term') matches anova(<list>)")
data(varespec, varechem)
perm <- shuffleSet(nrow(varespec), 49)
mod <- cca(varespec ~ Al + P + K, varechem)
anoterm <- anova(mod, by = "term", permutations = perm)
ano1df <- anova(mod, by = "onedf", permutations = perm)
expect_equal(permustats(anoterm)$permutations,
             permustats(ano1df)$permutations,
             info = "anova by 'term' and 'onedf' match")

## anova.ccanull:
## unconstrained model
mod <- ca(dune)
mat <- matrix(c(0, mod$CA$rank, 0, mod$tot.chi, 0, NA, NA, NA), nrow=2, ncol=4)
expect_equal(unname(as.matrix(anova(mod))), mat)
## overfitted model
expect_message(mod <- cca(dune ~ (.)^2, dune.env), "model is overfitted")
mat <- matrix(c(mod$CCA$rank, 0, mod$tot.chi, 0, NA, NA, NA, NA),
              nrow=2, ncol=4)
expect_equal(unname(as.matrix(anova(mod))), mat)

### cca as a weigthed RDA of Chi-square standardized data
if (packageVersion("vegan") >= "2.8") { # export ordConstrained
mod <- cca(dune ~ Moisture + Management, dune.env)
perm <- shuffleSet(nrow(dune), 49)
## see help(ordConstrained)
initWPCA <- function(Y, w) {
    Y <- as.matrix(Y)
    Y <- .Call("do_wcentre", Y, w)
    attr(Y, "RW") <- w
    attr(Y, "METHOD") <- "WPCA"
    Y
}
w <- rowSums(dune)/sum(dune)
modw <- ordConstrained(initWPCA(decostand(dune, "chi.sq"), w),
                       model.matrix(mod),
                       NULL, "pass")
class(modw) <- c("wrda", "cca")
ano <- anova(mod, permutations = perm)
anow <- anova(modw, permutations = perm)
expect_equal(permustats(ano)$permutations, permustats(anow)$permutations,
             info = "cca is weighted rda")
}

if (packageVersion("vegan") >= "2.8") { # parametric anova
### rda: parametric anova
data(varespec, varechem)
mod <- rda(varespec ~ Al + P + N, varechem)
mod0 <- lm(as.matrix(varespec) ~ Al + P + N, varechem)
ano <- anova(mod, by = "term", test = "parametric")
an0 <- anova(mod0, test="Spherical")
expect_equal(ano[["Pr(>F)"]], an0[["G-G Pr"]][-1],  # minus (Intercept)
             info = "parametric rda and mlm match")
## univariate model
mod <- rda(diversity(varespec) ~ Baresoil + Ca + P, varechem)
mod0 <- lm(diversity(varespec) ~ Baresoil + Ca + P, varechem)
expect_equal(anova(mod, test="parametric", by="term")$Pr, anova(mod0)$Pr,
             info = "univariate parametric rda and lm match")

### Parametric interface and forwarding through model-selection functions
interface.response <- dune.env$A1 + rep(c(-0.1, 0.1), length.out = nrow(dune))
interface.null <- rda(interface.response ~ 1, dune.env)
interface.full <- rda(interface.response ~ A1, dune.env)
interface.lm0 <- lm(interface.response ~ 1, dune.env)
interface.lm1 <- lm(interface.response ~ A1, dune.env)
expect_equal(anova(interface.full, test = "parametric")$Pr[1],
             anova(interface.lm1)$Pr[1])
expect_equal(anova(interface.null, interface.full, test = "parametric")$Pr[2],
             anova(interface.lm0, interface.lm1)$Pr[2])
expect_equal(anova(interface.full, test = "para"),
             anova(interface.full, test = "parametric"))
expect_equal(add1(interface.null, scope = formula(interface.full),
                  test = "parametric")$Pr,
             add1(interface.lm0, scope = formula(interface.lm1), test = "F")$Pr)
expect_equal(drop1(interface.full, test = "parametric")$Pr,
             drop1(interface.lm1, test = "F")$Pr)
interface.step <- ordistep(interface.null, scope = formula(interface.full),
                           direction = "forward", test = "parametric",
                           trace = FALSE)
expect_equal(attr(terms(interface.step), "term.labels"), "A1")
interface.r2step <- ordiR2step(interface.null, scope = formula(interface.full),
                               test = "parametric", trace = FALSE)
expect_equal(attr(terms(interface.r2step), "term.labels"), "A1")
expect_error(anova(interface.full, test = "F"))
expect_error(add1(interface.null, scope = formula(interface.full), test = "F"))
expect_error(drop1(interface.full, test = "F"))
expect_error(ordistep(interface.null, scope = formula(interface.full),
                      test = "F", trace = FALSE))
expect_error(ordiR2step(interface.null, scope = formula(interface.full),
                        test = "F", trace = FALSE))
}

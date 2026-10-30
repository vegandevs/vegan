### verification tests for anova.cca cases (incl. adonis2)

data(dune, dune.env)
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

### cca as a weigthed RDA of Chi-square standardized data
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

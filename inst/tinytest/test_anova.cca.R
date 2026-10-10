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
expect_silent(anova(mod0, by = "margin", test = "F"))

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
ano <- anova(mod, by = "term", test = "F")
an0 <- anova(mod0, test="Spherical")
expect_equal(ano[["Pr(>F)"]], an0[["G-G Pr"]][-1],  # minus (Intercept)
             info = "parametric rda and mlm match")
## univariate model
mod <- rda(diversity(varespec) ~ Baresoil + Ca + P, varechem)
mod0 <- lm(diversity(varespec) ~ Baresoil + Ca + P, varechem)
expect_equal(anova(mod, test="F", by="term")$Pr, anova(mod0)$Pr,
             info = "univariate parametric rda and lm match")

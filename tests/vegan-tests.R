### vegan-tests: unit tests for vegan functions

### This file contains unit tests for vegan functions. This file is
### run in R CMD check and its results are compared against previously
### saved results in vegan-tests.Rout.save. If you change tests, you
### must generate new vegan-tests.Rout.save in this directory.

### The current plan is that tests/ are not included in the CRAN
### release, but only in the development versin of vegan in R-Forge.

### The tests here are not intended for human reading. The tests need
### not be ecological or biologically meaningful, but they are only
### intended for testing strange arguments, protect against
### regressions and test correctness of results.

### The tests are in a single file and you should clean work space
### after the unit test. You should set random number seed (if needed)
### for comparison against vegan-tests.Rout.save, and you should
### remove the seed after the unit test. If possible, avoid very long
### lasting tests.

###<--- BEGIN TESTS --->
suppressPackageStartupMessages(require(vegan))
###<--- BEGIN anova.cca test --->
### anova.cca tests: should work with (1) subset, (2) missing values,
### (3) with functions of variables poly(A1,2), (4) variables in data
### frame attached or in data=, or (5) in working environment
set.seed(4711)
data(dune)
data(dune.env)
df <- dune.env
df$Management[c(1,5)] <- NA
## formula
fla <- as.formula("dune ~ Management + poly(A1, 2) + spno")
### variable in the .GlobalEnv
spno <- specnumber(dune)
### data= argument
## cca/rda
m <-  cca(fla, data=df,  na.action=na.exclude,  subset = Use != "Pasture" & spno > 7)
anova(m, perm=100)
anova(m, by="term", perm=100)
anova(m, by="margin", perm=100)
anova(m, by="axis", perm=100)
## capscale
p <- capscale(fla, data=df, na.action=na.exclude, subset = Use != "Pasture" & spno > 7)
anova(p, perm=100)
anova(p, by="term", perm=100)
anova(p, by="margin", perm=100)
anova(p, by="axis", perm=100)
## see that capscale can be updated and also works with 'dist' input
dis <- vegdist(dune)
p <- update(p, dis ~ .)
anova(p, perm=100)
anova(p, by="term", perm=100)
anova(p, by="margin", perm=100)
anova(p, by="axis", perm=100)
### attach()ed data frame instead of data=
attach(df)
q <- cca(fla, na.action = na.omit, subset = Use != "Pasture" & spno > 7)
anova(q, perm=100)
anova(q, by="term", perm=100)
anova(q, by="margin", perm=100)
anova(q, by="axis", perm=100)
### Check that constrained ordination functions can be embedded.
### The data.frame 'df' is still attach()ed.
foo <- function(bar, Y, X, ...)
{
    bar <- match.fun(bar)
    bar(Y ~ X, ...)
}
foo("cca", dune, Management, na.action = na.omit)
foo("rda", dune, Management, na.action = na.omit)
foo("capscale", dune, Management, dist="jaccard", na.action = na.omit)
foo("capscale", vegdist(dune), Management, na.action = na.omit)
### FIXME: foo("capscale", dune, Management, data=dune.env) fails!
###
detach(df)
### Check that statistics match in partial constrained ordination
m <- cca(dune ~ A1 + Moisture + Condition(Management), dune.env, subset = A1 > 3)
tab <- anova(m, by = "axis", perm.max = 100)
m
tab
all.equal(tab[,2], c(m$CCA$eig, m$CA$tot.chi), check.attributes=FALSE)
tab[nrow(tab),1] == m$CA$rank
## clean-up
rm(df, spno, fla, m, p, q, tab, dis, foo, .Random.seed)
### <--- END anova.cca test --->

### Sven Neulinger <sneulinger@ifam.uni-kiel.de> reported failures in
### partial analysis which (mostly) were fixed in r2087. Below his test.

set.seed(4711)
X <- matrix(rnorm(30*6), 30, 6)

A <- factor(rep(rep(c("a","b"), each=3),5))
B <- factor(rep(c("a","b","c"), 10))
## Sven Neulinger's tests used 'C' below, but that fails still now due
## to look-up order: function stats::C was found before matrix 'C'
CC <- factor(rep(c(1:5), each=6))

# partial db-RDA
cap.model.cond <- capscale(X ~ A + B + Condition(CC))
anova(cap.model.cond, by="axis", strata=CC)  # -> error pre r2287
anova(cap.model.cond, by="terms", strata=CC)  # -> error pre r2287

# db-RDA without conditional factor
cap.model <- capscale(X ~ A + B)
anova(cap.model, by="axis", strata=CC)  # -> no error
anova(cap.model, by="terms", strata=CC)  # -> no error

# partial RDA
rda.model.cond <- rda(X ~ A + B + Condition(CC))
anova(rda.model.cond, by="axis", strata=CC)  # -> no error
anova(rda.model.cond, by="terms", strata=CC)  # -> error pre r2287

# RDA without conditional factor
rda.model <- rda(X ~ A + B)
anova(rda.model, by="axis", strata=CC)  # -> no error
anova(rda.model, by="terms", strata=CC)  # -> no error
## clean.up
rm(X, A, B, CC, cap.model.cond, cap.model, rda.model.cond, rda.model)
### end Sven Neulinger's tests


### nestednodf: test case by Daniel Spitale in a comment to News on
### the release of vegan 1.17-6 in vegan.r-forge.r-project.org.
x <- c(1,0,1,1,1,1,1,1,0,0,0,1,1,1,0,1,1,0,0,0,1,1,0,0,0)
m1 <- matrix(x, nrow=5, ncol=5, byrow=FALSE)# as in Fig 2 Almeida-Neto et al 2008.
(nodf1 <- nestednodf(m1, order = FALSE, weighted = FALSE))
## Now the same matrix but with abundance data
x <- c(5,0,2,1,1,4,1,1,0,0,0,7,1,1,0,3,1,0,0,0,9,1,0,0,0)
m <- matrix(x, nrow=5, ncol=5, byrow=FALSE)
(nodfq <- nestednodf(m, order = FALSE, weighted = FALSE))
identical(nodf1, nodfq)
rm(x, m, m1, nodfq, nodf1)
### end nestednodf

### envfit & plot.envfit: latter failed if na.action resulted in only
### observation with a given factor level was removed. plot.envfit would
### fail with error about too long subscript
### fixed case where data presented to envfit also has extraneous levels
data(dune)
data(dune.env)
## add a new level to one of the factors
levels(dune.env$Management) <- c(levels(dune.env$Management), "foo")
## fit nMDS and envfit
set.seed(1)
mod <- metaMDS(dune)
ef <- envfit(mod, dune.env, permutations = 99)
plot(mod)
plot(ef, p.max = 0.1)
rm(mod, ef)
### end envfit & plot.envfit

### protest (& Procrustes analysis): Stability of the permutations and
### other results. 
data(mite)
mod <- rda(mite)
x <- scores(mod, display = "si", choices=1:6)
set.seed(4711)
xp <- x[sample(nrow(x)),]
pro <- protest(x, xp, permutations = 99)
pro
pro$t
rm(x, xp, pro)
### end protest

### Check that functions related to predict.rda work correctly for all
### constrained ordination methods.

### simulate.rda/cca/capscale: based on predict.* and the following
### should get back the data
data(dune, dune.env)
ind <- seq_len(nrow(dune))
target <- as.matrix(dune)
## rda
mod <- rda(dune ~ Condition(Moisture) + Management + A1, dune.env)
dat <- simulate(mod, indx = ind)
all.equal(dat, target, check.attributes = FALSE)
## cca
mod <- cca(dune ~ Condition(Moisture) + Management + A1, dune.env)
dat <- simulate(mod, indx = ind)
all.equal(dat, target, check.attributes = FALSE)
## capscale: Euclidean distances -- non-Euclidean distances have an
## imaginary component and will not give back the data.
d <- dist(dune)
mod <- capscale(d ~ Condition(Moisture) + Management + A1, dune.env)
dat <- simulate(mod, indx = ind)
all.equal(dat, d, check.attributes = FALSE)
## clean up
rm(ind, target, mod, dat, d)
### end simulate.*


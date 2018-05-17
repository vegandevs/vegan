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
anova(m, permutations=99)
anova(m, by="term", permutations=99) # failed before 2.5-0
anova(m, by="margin", permutations=99) # works since 2.5-0
anova(m, by="axis", permutations=99)
## capscale
p <- capscale(fla, data=df, na.action=na.exclude, subset = Use != "Pasture" & spno > 7)
anova(p, permutations=99)
anova(p, by="term", permutations=99) # failed before 2.5-0
anova(p, by="margin", permutations=99) # works since 2.5-0
anova(p, by="axis", permutations=99)
## see that capscale can be updated and also works with 'dist' input
dis <- vegdist(dune)
p <- update(p, dis ~ .)
anova(p, permutations=99)
anova(p, by="term", permutations=99) # failed before 2.5-0
anova(p, by="margin", permutations=99) # works since 2.5-0
anova(p, by="axis", permutations=99)
### attach()ed data frame instead of data=
attach(df)
q <- cca(fla, na.action = na.omit, subset = Use != "Pasture" & spno > 7)
anova(q, permutations=99)
## commented tests below fail in vegan 2.1-40 because number of
## observations changes
anova(q, by="term", permutations=99) # failed before 2.5-0
anova(q, by="margin", permutations=99) # works since 2.5-0
anova(q, by="axis", permutations=99)
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
foo("capscale", dune, Management, na.action = na.omit) ## fails in 2.2-1
###
detach(df)
### Check that statistics match in partial constrained ordination
m <- cca(dune ~ A1 + Moisture + Condition(Management), dune.env, subset = A1 > 3)
tab <- anova(m, by = "axis", permutations = 99)
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
## Sven Neulinger's tests failed still in 2.2-1, now due to look-up
## order: function stats::C was found before matrix 'C'. The test was
## OK when non-function name was used ('CC').
C <- factor(rep(c(1:5), each=6))

## partial db-RDA
cap.model.cond <- capscale(X ~ A + B + Condition(C))
anova(cap.model.cond, by="axis", strata=C)  # -> error pre r2287
anova(cap.model.cond, by="terms", strata=C)  # -> error pre r2287

## db-RDA without conditional factor
cap.model <- capscale(X ~ A + B)
anova(cap.model, by="axis", strata=C)  # -> no error
anova(cap.model, by="terms", strata=C)  # -> no error

# partial RDA
rda.model.cond <- rda(X ~ A + B + Condition(C))
anova(rda.model.cond, by="axis", strata=C)  # -> no error
anova(rda.model.cond, by="terms", strata=C)  # -> error pre r2287

# RDA without conditional factor
rda.model <- rda(X ~ A + B)
anova(rda.model, by="axis", strata=C)  # -> no error
anova(rda.model, by="terms", strata=C)  # -> no error
## clean.up
rm(X, A, B, C, cap.model.cond, cap.model, rda.model.cond, rda.model)
### end Sven Neulinger's tests

### Benedicte Bachelot informed us that several anova.cca* functions
### failed if community data name was the same as a function name: the
### function name was found first, and used instead ofa data. This
### seems to be related to the same problem that Sven Neulinger
### communicated, and his examples still faile if Condition or strata
### are function names. However, the following examples that failed
### should work now:

set.seed(4711)
cca <- dune
m <- cca(cca ~ ., dune.env)
anova(m, by="term")
m <- capscale(cca ~ ., dune.env)
anova(m, by="term")
rm(m, cca)

### end Benedicte Bachelot tests

### Richard Telford tweeted this example on 23/2/2015. Fails in 2.2-1,
### but should work now. Also issue #100 in github.com/vegandevs/vegan.
set.seed(4711)
data(dune, dune.env)
foo <- function(x, env) {
    m <- rda(x ~ Manure + A1, data = env)
    anova(m, by = "margin")
}
out <- lapply(dune, foo, env = dune.env)
out$Poatriv
rm(foo, out)
### end Richard Telford test


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
mod <- metaMDS(dune, trace = 0)
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
pro <- protest(x, xp, permutations = how(nperm = 99))
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

### test metaMDS works with long expression for comm
### originally reported to GLS by Richard Telford
data(varespec)
set.seed(1)
mod <- metaMDS(subset(varespec, select = colSums(varespec) > 0, subset = rowSums(varespec) > 0), trace=0)
mod
rm(mod)
### The above should run without error & last lines tests changes to the
### printed output

## dbrda tests

## the following three should be all equal
data(varespec, varechem)
(mr <- rda(varespec ~ Al + P + Condition(pH), varechem))
(md <- dbrda(varespec ~ Al + P + Condition(pH), varechem))
(mc <- capscale(varespec ~ Al + P + Condition(pH), varechem))
## the following two should be zero (within 1e-15)
p <- shuffleSet(nrow(varespec), 999)
all(abs(permustats(anova(mr, permutations=p))$permutations -
        permustats(anova(md, permutations=p))$permutations)
             < sqrt(.Machine$double.eps))

all(abs(permustats(anova(mr, permutations=p))$permutations -
        permustats(anova(mc, permutations=p))$permutations)
             < sqrt(.Machine$double.eps))
## eigenvals returns a list now (>= 2.5-0)
data(varespec, varechem)
mod <- cca(varespec ~ Al + P + Condition(pH), varechem)
ev <- summary(eigenvals(mod))
stopifnot(inherits(ev, "matrix"))
stopifnot(!is.list(ev))
ev

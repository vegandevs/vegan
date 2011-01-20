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
require(vegan, quiet = TRUE)
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
### attach()ed data frame instead of data=
attach(df)
q <- cca(fla, na.action = na.omit, subset = Use != "Pasture" & spno > 7)
anova(q, perm=100)
anova(q, by="term", perm=100)
anova(q, by="margin", perm=100)
anova(q, by="axis", perm=100)
detach(df)
## clean-up
rm(df, spno, fla, m, p, q, .Random.seed)
### <--- END anova.cca test --->

## test the stability of cca object and its support functions
suppressPackageStartupMessages(require(vegan))
set.seed(4711)
## models
data(dune, dune.env)
mcca <- cca(dune ~ Condition(Management) + Manure + A1, dune.env)
mrda <- rda(dune ~ Condition(Management) + Manure + A1, dune.env)
mcap <- capscale(dune ~ Condition(Management) + Manure + A1, dune.env,
          dist = "bray")
mdb <- dbrda(dune ~ Condition(Management) + Manure + A1, dune.env,
             dist = "bray")
## general appearance
mcca
mrda
mcap
mdb
## names
sort(names(mcca))
sort(names(mrda))
sort(names(mcap))
sort(names(mdb))
## overall results
str(mcca, max.level = 2)
str(mrda, max.level = 2)
str(mcap, max.level = 2)
str(mdb, max.level=2)

## diagnostics

head(goodness(mcca, display = "sites"))
head(goodness(mrda, display = "sites"))
head(goodness(mcap, display = "sites"))
head(goodness(mdb, display="sites"))

head(inertcomp(mcca))
head(inertcomp(mrda))
head(inertcomp(mcap, display="sites"))
head(inertcomp(mdb, display = "sites"))

intersetcor(mcca)
intersetcor(mrda)
intersetcor(mcap)
intersetcor(mdb)

tolerance(mcca)

vif.cca(mcca)
vif.cca(mrda)
vif.cca(mcap)
vif.cca(mdb)

alias(mcca)
alias(mrda)
alias(mcap)
alias(mdb)

## basic statistic

coef(mcca)
coef(mrda)
coef(mcap)
coef(mdb)

eigenvals(mcca)
eigenvals(mrda)
eigenvals(mcap)
eigenvals(mdb)

nobs(mcca)
nobs(mrda)
nobs(mcap)
nobs(mdb)

RsquareAdj(mcca)
RsquareAdj(mrda)
RsquareAdj(mcap)
RsquareAdj(mdb)

## casting

as.mlm(mcca)
as.mlm(mrda)
as.mlm(mcap)
as.mlm(mdb)

head(model.frame(mcca))
head(model.frame(mrda))
head(model.frame(mcap))
head(model.frame(mdb))

## testing and model building -- NB, most anova.cca tests and model
## building are in are in vegan-tests.R and in Examples/vegan-Ex.R

deviance(mcca)
deviance(mrda)
deviance(mcap)
deviance(mdb)

per <- shuffleSet(nrow(dune), 49)
permutest(mcca, per)
permutest(mrda, per)
permutest(mcap, per)
permutest(mdb, per)

drop1(mcca)
drop1(mrda)
drop1(mcap)
drop1(mdb)

## the following do not work with partial models

mcca <- cca(dune ~ Management + Manure + A1, dune.env)
mrda <- rda(dune ~ Management + Manure + A1, dune.env)
mcap <- capscale(dune ~ Management + Manure + A1, dune.env, dist = "bray")
mdb <- dbrda(dune ~ Management + Manure + A1, dune.env, dist = "bray")

head(calibrate(mcca))
head(calibrate(mrda))
head(calibrate(mcap))
head(calibrate(mdb))

head(calibrate(mcca, newdata=dune[11:15,]))
head(calibrate(mrda, newdata=dune[11:15,]))




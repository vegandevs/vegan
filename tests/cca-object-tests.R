## test the stability of cca object and its support functions
suppressPackageStartupMessages(require(vegan))
set.seed(4711)
## models
data(dune, dune.env)
mcca <- cca(dune ~ Condition(Management) + Manure + A1, dune.env)
mrda <- rda(dune ~ Condition(Management) + Manure + A1, dune.env)
mrda1 <- rda(dune ~ Condition(Management) + Manure + A1, dune.env, scale=TRUE)
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
head(goodness(mrda1, display = "sites"))
head(goodness(mcap, display = "sites"))
head(goodness(mdb, display="sites"))

head(inertcomp(mcca))
head(inertcomp(mrda))
head(inertcomp(mrda1))
head(inertcomp(mcap, display="sites"))
head(inertcomp(mdb, display = "sites"))

intersetcor(mcca)
intersetcor(mrda)
intersetcor(mrda1)
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
coef(mrda1)
coef(mcap)
coef(mdb)

eigenvals(mcca)
eigenvals(mrda)
eigenvals(mrda1)
eigenvals(mcap)
eigenvals(mdb)

nobs(mcca)
nobs(mrda)
nobs(mcap)
nobs(mdb)

RsquareAdj(mcca)
RsquareAdj(mrda)
RsquareAdj(mrda1)
RsquareAdj(mcap)
RsquareAdj(mdb)

## casting

as.mlm(mcca)
as.mlm(mrda)
as.mlm(mrda1)
as.mlm(mcap)
as.mlm(mdb)

head(model.frame(mcca))
head(model.frame(mrda))
head(model.frame(mrda1))
head(model.frame(mcap))
head(model.frame(mdb))

## testing and model building -

deviance(mcca)
deviance(mrda)
deviance(mrda1)
deviance(mcap)
deviance(mdb)

per <- shuffleSet(nrow(dune), 49)
permutest(mcca, per)
permutest(mrda, per)
permutest(mrda1, per)
permutest(mcap, per)
permutest(mdb, per)

drop1(mcca, test="permutation", permutations=per)
drop1(mrda, test="permutation", permutations=per)
drop1(mrda1, test="permutation", permutations=per)
drop1(mcap, test="permutation", permutations=per)
drop1(mdb, test="permutation", permutations=per)

anova(mcca, permutations = per)
anova(mrda, permutations = per)
anova(mrda1, permutations = per)
anova(mcap, permutations = per)
anova(mdb, permutations = per)

anova(mcca, permutations = per, by="term")
anova(mrda, permutations = per, by="term")
anova(mrda1, permutations = per, by="term")
anova(mcap, permutations = per, by="term")
anova(mdb, permutations = per, by="term")

anova(mcca, permutations = per, by="margin")
anova(mrda, permutations = per, by="margin")
anova(mrda1, permutations = per, by="margin")
anova(mcap, permutations = per, by="margin")
anova(mdb, permutations = per, by="margin")

anova(mcca, permutations = per, by="axis")
anova(mrda, permutations = per, by="axis")
anova(mrda1, permutations = per, by="axis")
anova(mcap, permutations = per, by="axis")
anova(mdb, permutations = per, by="axis")

## the following do not work with partial models

mcca <- cca(dune ~ Management + Manure + A1, dune.env)
mrda <- rda(dune ~ Management + Manure + A1, dune.env)
mrda1 <- rda(dune ~ Management + Manure + A1, dune.env, scale = TRUE)
mcap <- capscale(dune ~ Management + Manure + A1, dune.env, dist = "bray")
mdb <- dbrda(dune ~ Management + Manure + A1, dune.env, dist = "bray")

head(calibrate(mcca))
head(calibrate(mrda))
head(calibrate(mrda1))
head(calibrate(mcap))
head(calibrate(mdb))

head(calibrate(mcca, newdata=dune[11:15,]))
head(calibrate(mrda, newdata=dune[11:15,]))
head(calibrate(mrda1, newdata=dune[11:15,]))




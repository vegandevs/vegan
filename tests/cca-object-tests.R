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
mancap <- capscale(dune ~ Condition(Management) + Manure + A1, dune.env,
          dist = "manhattan")
mandb <- dbrda(dune ~ Condition(Management) + Manure + A1, dune.env,
               dist = "manhattan")
## 0-rank constraints
m0cca <- cca(dune ~ Condition(Management) + Management, dune.env)

## general appearance
mcca
mrda
mcap
mdb
mancap
mandb
m0cca
## names
sort(names(mcca))
sort(names(mrda))
sort(names(mcap))
sort(names(mdb))
## overall results
head(summary(mcca))
head(summary(mrda))
head(summary(mrda1))
head(summary(mcap))
head(summary(mdb))
head(summary(mancap))
head(summary(mandb))
head(summary(m0cca))
## diagnostics
hatvalues(mcca)
hatvalues(mrda)
hatvalues(mandb)

head(cooks.distance(mcca))
head(cooks.distance(mrda))
head(cooks.distance(mrda1))
head(cooks.distance(mcap, "canoco"))
head(cooks.distance(mdb, "canoco"))
head(cooks.distance(mancap, "canoco"))
head(cooks.distance(mandb, "canoco"))
head(cooks.distance(m0cca))

head(goodness(mcca, display = "sites"))
head(goodness(mrda, display = "sites"))
head(goodness(mrda1, display = "sites"))
## head(goodness(mcap, display = "sites")) # currently disabled
## head(goodness(mdb, display="sites"))  # not implemented for partial dbrda
## head(goodness(mancap, display="sites")) # currently disabled
## head(goodness(mandb, display="sites")) # not implemneted for partial dbrda
head(goodness(m0cca))

head(inertcomp(mcca))
head(inertcomp(mrda))
head(inertcomp(mrda1))
head(inertcomp(mcap, display="sites"))
head(inertcomp(mdb, display = "sites"))
head(inertcomp(mancap, display = "sites"))
head(inertcomp(mandb, display = "sites"))
head(inertcomp(m0cca))

zapsmall(intersetcor(mcca))
zapsmall(intersetcor(mrda))
zapsmall(intersetcor(mrda1))
zapsmall(intersetcor(mcap))
zapsmall(intersetcor(mdb))
zapsmall(intersetcor(mancap))
zapsmall(intersetcor(mandb))

tolerance(mcca)
tolerance(m0cca)

vif.cca(mcca)
vif.cca(mrda)
vif.cca(mcap)
vif.cca(mdb)
vif.cca(m0cca)

alias(mcca)
alias(mrda)
alias(mcap)
alias(mdb)
alias(m0cca)

## basic statistic

coef(mcca)
coef(mrda)
coef(mrda1)
coef(mcap)
coef(mdb)
coef(m0cca)

eigenvals(mcca)
eigenvals(mrda)
eigenvals(mrda1)
eigenvals(mcap)
eigenvals(mdb)
eigenvals(mancap)
eigenvals(mandb)
eigenvals(m0cca)
eigenvals(m0cca, model = "constrained")

nobs(mcca)
nobs(mrda)
nobs(mcap)
nobs(mdb)
nobs(m0cca)

RsquareAdj(mcca)
RsquareAdj(mrda)
RsquareAdj(mrda1)
RsquareAdj(mcap)
RsquareAdj(mdb)
RsquareAdj(m0cca)

head(model.frame(mcca))
head(model.frame(mrda))
head(model.frame(mrda1))
head(model.frame(mcap))
head(model.frame(mdb))
head(model.frame(m0cca))

## testing and model building -

deviance(mcca)
deviance(mrda)
deviance(mrda1)
deviance(mcap)
deviance(mdb)
deviance(m0cca)

per <- shuffleSet(nrow(dune), 49)
permutest(mcca, per)
permutest(mrda, per)
permutest(mrda1, per)
permutest(mcap, per)
permutest(mdb, per)
permutest(mancap, per)
permutest(mandb, per)
permutest(m0cca, per)

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
anova(mancap, permutations = per)
anova(mandb, permutations = per)
anova(m0cca, permutations = per)

anova(mcca, permutations = per, by="term")
anova(mrda, permutations = per, by="term")
anova(mrda1, permutations = per, by="term")
anova(mcap, permutations = per, by="term")
anova(mdb, permutations = per, by="term")
anova(mancap, permutations = per, by="term")
anova(mandb, permutations = per, by="term")

anova(mcca, permutations = per, by="margin")
anova(mrda, permutations = per, by="margin")
anova(mrda1, permutations = per, by="margin")
anova(mcap, permutations = per, by="margin")
anova(mdb, permutations = per, by="margin")
anova(mancap, permutations = per, by="margin")
anova(mandb, permutations = per, by="margin")

anova(mcca, permutations = per, by="axis")
anova(mrda, permutations = per, by="axis")
anova(mrda1, permutations = per, by="axis")
anova(mcap, permutations = per, by="axis")
anova(mdb, permutations = per, by="axis")
anova(mancap, permutations = per, by="axis")
anova(mandb, permutations = per, by="axis")

## the following do not all work with partial models

mcca <- cca(dune ~ Management + Manure + A1, dune.env)
mrda <- rda(dune ~ Management + Manure + A1, dune.env)
mrda1 <- rda(dune ~ Management + Manure + A1, dune.env, scale = TRUE)
mcap <- capscale(dune ~ Management + Manure + A1, dune.env, dist = "bray")
mdb <- dbrda(dune ~ Management + Manure + A1, dune.env, dist = "bray")
mancap <- capscale(dune ~ Management + Manure + A1, dune.env, dist = "man")
mandb <- dbrda(dune ~ Management + Manure + A1, dune.env, dist = "man")

head(calibrate(mcca))
head(calibrate(mrda))
head(calibrate(mrda1))
head(calibrate(mcap))
head(calibrate(mdb))
head(calibrate(mancap))
head(calibrate(mandb))

head(calibrate(mcca, newdata=dune[11:15,]))
head(calibrate(mrda, newdata=dune[11:15,]))
head(calibrate(mrda1, newdata=dune[11:15,]))

head(predict(mcca, newdata = dune.env))
predict(mrda, newdata = dune.env[1:4,])
predict(mrda1, newdata = dune.env[1:4,])
predict(mcap, newdata = dune.env[1:4,])
predict(mdb, newdata = dune.env[1:4,])
predict(mancap, newdata = dune.env[1:4,])
predict(mandb, newdata = dune.env[1:4,])

predict(mcca, newdata = dune[1:4,], type="wa")
predict(mrda, newdata = dune[1:4,], type="wa")
predict(mrda1, newdata = dune[1:4,], type="wa")

predict(mcca, newdata = dune[,1:4], type="sp")
predict(mrda, newdata = dune[,1:4], type="sp")
predict(mrda1, newdata = dune[,1:4], type="sp")

predict(mcca, newdata = dune.env[1:4,], type="lc")
predict(mrda, newdata = dune.env[1:4,], type="lc")
predict(mrda1, newdata = dune.env[1:4,], type="lc")
predict(mcap, newdata = dune.env[1:4,], type="lc")
predict(mdb, newdata = dune.env[1:4,], type="lc")
predict(mancap, newdata = dune.env[1:4,], type="lc")
predict(mandb, newdata = dune.env[1:4,])


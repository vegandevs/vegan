## test the stability of cca object and its support functions
suppressPackageStartupMessages(require(vegan))
suppressPackageStartupMessages(require(parallel))
set.seed(4711)
op <- options(digits=5)
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
## univariate model, rank-1 constraint
H <- diversity(dune)
m1rda <- rda(H ~  A1, dune.env)

## general appearance
mcca
mrda
mcap
mdb
mancap
mandb
m0cca
m1rda
## names
sort(names(mcca))
sort(names(mrda))
sort(names(mcap))
sort(names(mdb))

## diagnostics
hatvalues(mcca)
hatvalues(mrda)
hatvalues(mandb)
hatvalues(m1rda)

zapsmall(head(cooks.distance(mcca)))
zapsmall(head(cooks.distance(mrda)))
zapsmall(head(cooks.distance(mrda1)))
zapsmall(head(cooks.distance(mcap, "canoco")))
zapsmall(head(cooks.distance(mdb, "canoco")))
zapsmall(head(cooks.distance(mancap, "canoco")))
zapsmall(head(cooks.distance(mandb, "canoco")))
zapsmall(head(cooks.distance(m1rda)))

head(goodness(mcca, display = "sites"))
head(goodness(mrda, display = "sites"))
head(goodness(mrda1, display = "sites"))
head(goodness(m1rda, display = "sites"))
goodness(m1rda)       # fails in 2.5-3
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
zapsmall(head(inertcomp(m0cca)))
inertcomp(m1rda)

abs(zapsmall(intersetcor(mcca)))
abs(zapsmall(intersetcor(mrda)))
abs(zapsmall(intersetcor(mrda1)))
abs(zapsmall(intersetcor(mcap)))
abs(zapsmall(intersetcor(mdb)))
abs(zapsmall(intersetcor(mancap)))
abs(zapsmall(intersetcor(mandb)))
abs(zapsmall(intersetcor(m1rda)))

tolerance(mcca)
tolerance(m0cca)

vif.cca(mcca)
vif.cca(mrda)
vif.cca(mcap)
vif.cca(mdb)

alias(mcca)
alias(mrda)
alias(mcap)
alias(mdb)

## basic statistic

abs(coef(mcca))
abs(coef(mrda))
abs(coef(mrda1))
abs(coef(mcap))
abs(coef(mdb))
abs(coef(m1rda))

eigenvals(mcca)
eigenvals(mrda)
eigenvals(mrda1)
eigenvals(mcap)
eigenvals(mdb)
eigenvals(mancap)
eigenvals(mandb)
eigenvals(m0cca)
eigenvals(m0cca, model = "constrained")
eigenvals(m1rda)

nobs(mcca)
nobs(mrda)
nobs(mcap)
nobs(mdb)
nobs(m0cca)
nobs(m1rda)

RsquareAdj(mcca)
RsquareAdj(mrda)
RsquareAdj(mrda1)
RsquareAdj(mcap)
RsquareAdj(mdb)
RsquareAdj(m1rda)

head(model.frame(mcca))
head(model.frame(mrda))
head(model.frame(mrda1))
head(model.frame(mcap))
head(model.frame(mdb))
head(model.frame(m0cca))
head(model.frame(m1rda))

## testing and model building -

deviance(mcca)
deviance(mrda)
deviance(mrda1)
deviance(mcap)
deviance(mdb)
deviance(m0cca)
deviance(m1rda)

per <- shuffleSet(nrow(dune), 49)
permutest(mcca, per)
permutest(mrda, per)
permutest(mrda1, per)
permutest(mcap, per)
permutest(mdb, per)
permutest(mancap, per)
permutest(mandb, per)
permutest(m1rda, per)

drop1(mcca,  test = "permutation", permutations = per)
drop1(mrda,  test = "permutation", permutations = per)
drop1(mrda1, test = "permutation", permutations = per)
drop1(mcap,  test = "permutation", permutations = per)
drop1(mdb,   test = "permutation", permutations = per)
drop1(m1rda, test = "permutation", permutations = per)

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
anova(m1rda, permutations = per, by="term")

anova(mcca, permutations = per, by="margin")
anova(mrda, permutations = per, by="margin")
anova(mrda1, permutations = per, by="margin")
anova(mcap, permutations = per, by="margin")
anova(mdb, permutations = per, by="margin")
anova(mancap, permutations = per, by="margin")
anova(mandb, permutations = per, by="margin")
anova(m1rda, permutations = per, by="margin")

anova(mcca, permutations = per, by="axis")
anova(mrda, permutations = per, by="axis")
anova(mrda1, permutations = per, by="axis")
anova(mcap, permutations = per, by="axis")
anova(mdb, permutations = per, by="axis")
anova(mancap, permutations = per, by="axis")
anova(mandb, permutations = per, by="axis")
anova(m1rda, permutations = per, by="axis")

## permutation tests in parallel
clust <- makeCluster(2) # socket cluster: the only one that works in Windows

permutest(mcca,   per, parallel = clust) # use socket cluster
permutest(mrda,   per, parallel = clust) # use socket cluster
permutest(mrda1,  per, parallel = clust) # use socket cluster
permutest(mcap,   per, parallel = clust) # use socket cluster
permutest(mdb,    per, parallel = clust) # use socket cluster
permutest(mancap, per, parallel = clust) # use socket cluster
permutest(mandb,  per, parallel = clust) # use socket cluster
permutest(m1rda,  per, parallel = clust) # use socket cluster


drop1(mcca,  test = "permutation", permutations = per, parallel = clust)
drop1(mrda,  test = "permutation", permutations = per, parallel = clust)
drop1(mrda1, test = "permutation", permutations = per, parallel = clust)
drop1(mcap,  test = "permutation", permutations = per, parallel = clust)
drop1(mdb,   test = "permutation", permutations = per, parallel = clust)
drop1(m1rda, test = "permutation", permutations = per, parallel = clust)

anova(mcca,   permutations = per, parallel = clust)
anova(mrda,   permutations = per, parallel = clust)
anova(mrda1,  permutations = per, parallel = clust)
anova(mcap,   permutations = per, parallel = clust)
anova(mdb,    permutations = per, parallel = clust)
anova(mancap, permutations = per, parallel = clust)
anova(mandb,  permutations = per, parallel = clust)
anova(m0cca,  permutations = per, parallel = clust)

anova(mcca,   permutations = per, by = "term", parallel = clust)
anova(mrda,   permutations = per, by = "term", parallel = clust)
anova(mrda1,  permutations = per, by = "term", parallel = clust)
anova(mcap,   permutations = per, by = "term", parallel = clust)
anova(mdb,    permutations = per, by = "term", parallel = clust)
anova(mancap, permutations = per, by = "term", parallel = clust)
anova(mandb,  permutations = per, by = "term", parallel = clust)
anova(m1rda,  permutations = per, by = "term", parallel = clust)

anova(mcca,   permutations = per, by = "margin", parallel = clust)
anova(mrda,   permutations = per, by = "margin", parallel = clust)
anova(mrda1,  permutations = per, by = "margin", parallel = clust)
anova(mcap,   permutations = per, by = "margin", parallel = clust)
anova(mdb,    permutations = per, by = "margin", parallel = clust)
anova(mancap, permutations = per, by = "margin", parallel = clust)
anova(mandb,  permutations = per, by = "margin", parallel = clust)
anova(m1rda,  permutations = per, by = "margin", parallel = clust)

anova(mcca,   permutations = per, by = "axis", parallel = clust)
anova(mrda,   permutations = per, by = "axis", parallel = clust)
anova(mrda1,  permutations = per, by = "axis", parallel = clust)
anova(mcap,   permutations = per, by = "axis", parallel = clust)
anova(mdb,    permutations = per, by = "axis", parallel = clust)
anova(mancap, permutations = per, by = "axis", parallel = clust)
anova(mandb,  permutations = per, by = "axis", parallel = clust)
anova(m1rda,  permutations = per, by = "axis", parallel = clust)

# stop the cluster as we are finished
stopCluster(clust)

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
## head(calibrate(m1rda, newdata=dune[11:15,]))## fails

head(predict(mcca, newdata = dune.env))
predict(mrda, newdata = dune.env[1:4,])
predict(mrda1, newdata = dune.env[1:4,])
predict(mcap, newdata = dune.env[1:4,])
predict(mdb, newdata = dune.env[1:4,])
predict(mancap, newdata = dune.env[1:4,])
predict(mandb, newdata = dune.env[1:4,])
predict(m1rda, newdata = dune.env[1:4,])

## the sign is arbitrary
abs(predict(mcca, newdata = dune[1:4,], type="wa"))
abs(predict(mrda, newdata = dune[1:4,], type="wa"))
abs(predict(mrda1, newdata = dune[1:4,], type="wa"))

## the sign is arbitrary
abs(predict(mcca, newdata = dune[,1:4], type="sp"))
abs(predict(mrda, newdata = dune[,1:4], type="sp"))
abs(predict(mrda1, newdata = dune[,1:4], type="sp"))

abs(predict(mcca, newdata = dune.env[1:4,], type="lc"))
abs(predict(mrda, newdata = dune.env[1:4,], type="lc"))
abs(predict(mrda1, newdata = dune.env[1:4,], type="lc"))
abs(predict(mcap, newdata = dune.env[1:4,], type="lc"))
abs(predict(mdb, newdata = dune.env[1:4,], type="lc"))
abs(predict(mancap, newdata = dune.env[1:4,], type="lc"))
abs(predict(m1rda, newdata = dune.env[1:4,], type="lc"))
abs(predict(mandb, newdata = dune.env[1:4,]))
## reset
options(op)

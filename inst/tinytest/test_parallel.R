### tests for parallel processing: all choices should work and
### replicate same result sequences

### parallel=1 will always be non-parallel (sequential), and parallel
### = SOCK will use pre-defined socket cluster, but actions with
### parallel = <n> will depend on OS. In unix-alikes (incl. macOS),
### the command will fork processes, but in Windows a socket cluster
### is created in the function and closed on.exit. Parallel code was
### added in vegan 2.2-0 (2014). I had Linux desktop and Mac laptop,
### and fork clusters worked in both, but did not work in
### Windows. Support for socket clusters was added for Windows users
### with scanty testing. In 2014 fork clusters were clearly faster
### than sockets (which had larger overhead of starting the
### cluster). Now socket clusters are faster in ARM-based Mac, and
### this decision may be reconsidered.

### setup: CRAN only allows 2 parallel processes in their tests
NCPU <- min(parallel::detectCores(), 2)
SOCK <- parallel::makeCluster(NCPU)

### anosim: permutations can run in parallel
data(dune, dune.env)
perm <- shuffleSet(nrow(dune), 24)
p1 <- anosim(dune, dune.env$Management, permutations = perm, parallel = 1)
pfork <- anosim(dune, dune.env$Management, permutations = perm,
                parallel = NCPU)
psock <- anosim(dune, dune.env$Management, permutations = perm,
                parallel = SOCK)
expect_equivalent(pfork$perm, p1$perm)
## R 4.1 drops names in socket cluster in git tests
expect_equivalent(psock$perm, p1$perm)

## bioenv: models can be compared in parallel
data(varespec, varechem)
p1 <- bioenv(varespec ~ Al + P + K + N + Humdepth + pH, varechem,
             parallel = 1)
pfork <- bioenv(varespec ~ Al + P + K + N + Humdepth + pH, varechem,
             parallel = NCPU)
psock <- bioenv(varespec ~ Al + P + K + N + Humdepth + pH, varechem,
                parallel = SOCK)
expect_equal(pfork$models, p1$models)
expect_equal(psock$models, p1$models)

## cascadeKM: cannot be replicated with parallel processing since
## stats::kmeans() starts with unpredictable RNG seeds in parallel
## processes. The analysis cannot be replicated with parallel
## processes (same parallel processes with same seeds
## differ). Non-parallel cascadeKM can be replicated.
## Issue #771: socket clusters failed before 2.7-5.
data(dune)
set.seed(4711); p0 <- cascadeKM(dune, 2, 6, iter = 10, parallel = 1)
set.seed(4711); p1 <- cascadeKM(dune, 2, 6, iter = 10, parallel = 1)
expect_equal(p0, p1)
## WARNING: parallel runs cannot be replicated!
set.seed(4711); ppar0 <- cascadeKM(dune, 2, 6, iter = 10, parallel = SOCK)
set.seed(4711); ppar1 <- cascadeKM(dune, 2, 6, iter = 10, parallel = SOCK)
expect_false(isTRUE(all.equal(ppar0, ppar1)))

## estaccumR
data(BCI)
perm <- shuffleSet(nrow(BCI), 20)
p1 <- estaccumR(BCI, permutations = perm, parallel = 1)
pfork <- estaccumR(BCI, permutations = perm, parallel = NCPU)
psock <- estaccumR(BCI, permutations = perm, parallel = SOCK)
expect_equal(pfork, p1)
expect_equal(psock, p1)

## mantel & mantel.partial
data(mite, mite.env, mite.xy)
perm <- shuffleSet(nrow(mite), 24)
d <- vegdist(mite)
denv <- vegdist(mite.env[,1:2], "mahalanobis")
dgeo <- dist(mite.xy)
p1 <- mantel(d, denv, permutations = perm, parallel = 1)
pfork <- mantel(d, denv, permutations = perm, parallel = NCPU)
psock <- mantel(d, denv, permutations = perm, parallel = SOCK)
expect_equal(pfork$perm, p1$perm)
expect_equal(psock$perm, p1$perm)
p1 <- mantel.partial(d, denv, dgeo, permutations = perm, parallel = 1)
pfork <- mantel.partial(d, denv, dgeo, permutations = perm, parallel = NCPU)
psock <- mantel.partial(d, denv, dgeo, permutations = perm, parallel = SOCK)
expect_equal(pfork$perm, p1$perm)
expect_equal(psock$perm, p1$perm)

## metaMDS: metaMDSiter can run NMDS iterations in blocks of size 'parallel'
data(dune)
d <- vegdist(dune)
set.seed(4711); p1 <- metaMDS(d, trace = 0, trymax = 10, parallel = 1)
set.seed(4711); pfork <- metaMDS(d, trace = 0, trymax = 10, parallel = NCPU)
set.seed(4711); psock <- metaMDS(d, trace = 0, trymax = 10, parallel = SOCK)
items <- c("points", "stress", "bestry", "species")
expect_equal(pfork[items], p1[items])
expect_equal(psock[items], p1[items])

## mrpp
data(dune, dune.env)
perm <- shuffleSet(nrow(dune), 24)
p1 <- mrpp(dune, dune.env$Management, permutations = perm, parallel = 1)
pfork <- mrpp(dune, dune.env$Management, permutations = perm, parallel = NCPU)
psock <- mrpp(dune, dune.env$Management, permutations = perm, parallel = SOCK)
expect_equal(pfork$boot.deltas, p1$boot.deltas)
expect_equal(psock$boot.deltas, p1$boot.deltas)

## oecosimu can evaluate nestfun in parallel (this is pretty useless
## unless nestfun is very expensive). The test below actually checks
## null models are replicated; it would be possible to give nullmodel
## object as first argument to only test evaluation of nestfun.
data(sipoo)
set.seed(4711); p1 <- oecosimu(sipoo, nestedchecker, "r0", parallel = 1)
set.seed(4711); pfork <- oecosimu(sipoo, nestedchecker, "r0", parallel = NCPU)
set.seed(4711); psock <- oecosimu(sipoo, nestedchecker, "r0", parallel = SOCK)
expect_equivalent(pfork, p1)
expect_equivalent(psock, p1)

## ordiareatest (this had glitches pre-2.7-5: issue #774)
data(dune, dune.env)
perm <- shuffleSet(nrow(dune), 24)
ord <- cca(dune)
p1 <- ordiareatest(ord, dune.env$Management, permutations = perm, parallel = 1)
pfork <- ordiareatest(ord, dune.env$Management, permutations = perm,
                      parallel = NCPU)
psock <- ordiareatest(ord, dune.env$Management, permutations = perm,
                      parallel = SOCK)
expect_equal(pfork$permutations, p1$permutations)
expect_equal(psock$permutations, p1$permutations)

## permutest.betadisper. Failed pre-2.5.7 with socket clusters: issue #369.
data(dune, dune.env)
perm <- shuffleSet(nrow(dune), 24)
mod <- betadisper(vegdist(dune), dune.env$Management)
p1 <- permutest(mod, permutations = perm, pairwise = TRUE, parallel = 1)
pfork <- permutest(mod, permutations = perm, pairwise = TRUE, parallel = NCPU)
psock <- permutest(mod, permutations = perm, pairwise = TRUE, parallel = SOCK)
expect_equal(pfork, p1)
expect_equal(psock, p1)

## permutest.cca: anova.cca, drop1.cca, add1.cca, Rsquare.adj.cca
## delegate this method. Failed in socket clusters pre-2.5-2: issue
## #276.
data(dune, dune.env)
perm <- shuffleSet(nrow(dune), 24)
mod <- cca(dune ~ Management, dune.env)
p1 <- permutest(mod, by = "onedf", permutations = perm, parallel = 1)
pfork <- permutest(mod, by = "onedf", permutations = perm, parallel = NCPU)
psock <- permutest(mod, by = "onedf", permutations = perm, parallel = SOCK)
expect_equal(pfork$F.perm, p1$F.perm)
expect_equal(psock$F.perm, p1$F.perm)

### close socket cluster
parallel::stopCluster(SOCK)


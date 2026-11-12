### tests for running and handling NMDS basically metaMDS -> monoMDS

data(varespec, varechem)
expect_silent(m <- metaMDS(varespec, trymax=1, trace=3))
expect_silent(d <- metaMDSredist(m))
## monoMDS options
init <- wcmdscale(d) # starting position
expect_silent(mo <- monoMDS(d, init, model = "local", weakties = FALSE))
expect_silent(mlin <- monoMDS(d, init, model = "linear", stress = 2))
expect_silent(mhyb <- monoMDS(d, init, model = "hybrid", scaling = FALSE, pc = FALSE))
## goodness
expect_equal(sqrt(sum(goodness(m)^2)), m$stress)
expect_equal(sqrt(sum(goodness(mo)^2)), mo$stress)
expect_equal(sqrt(sum(goodness(mlin)^2)), mlin$stress)
expect_equal(sqrt(sum(goodness(mhyb)^2)), mhyb$stress)
expect_silent(m2 <- update(m, stress=2, trace=0)) # stress 2
expect_equal(sqrt(sum(goodness(m2)^2)), m2$stress)
expect_silent(mo2 <- update(mo, stress=2))
expect_equal(sqrt(sum(goodness(mo2)^2)), mo2$stress)
## stressplot: check data
stress1 <- function(x) sqrt(sum((x$y-x$yf)^2)/sum(x$y^2))
expect_silent(st <- stressplot(m))
expect_equal(stress1(st), m$stress)
expect_silent(st <- stressplot(mhyb))
expect_equal(stress1(st), mhyb$stress)
## add species scores to monoMDS
expect_length(scores(mo, "species"), 0) # no species scores
sppscores(mo) <- varespec
expect_equal(dim(scores(mo, "species")), c(ncol(varespec), mo$ndim))
## check calculation of species scores & adding them again
m0 <- m
m0$species <- NULL
expect_length(scores(m0, "species"), 0) # no species scores
sppscores(m0) <- eval(parse(text = m0$data)) # data as in metaMDS
expect_equal(scores(m0, "species"), scores(m, "species"))

if (packageVersion("vegan") >= "2.8") { # engine
## check engine and stepacross
data(dune) # varespec has no sites that 'noshare' species
expect_silent(miso <- metaMDS(dune, trymax = 0, trace = 0, noshare = TRUE,
                              engine = MASS::isoMDS))
expect_equal(sqrt(sum(goodness(miso)^2)), miso$stress)
## isoMDS uses 100*stress1
expect_equal(stress1(stressplot(miso)) * 100, miso$stress)
}

## Rotate to external variable
ef0 <- envfit(m ~ Al + Humdepth + Ca, varechem)
expect_silent(mAl <- MDSrotate(m, varechem$Al))
ef <- envfit(mAl ~ Al + Humdepth + Ca, varechem)
expect_equal(ef$vectors$r, ef0$vectors$r) # unchanged R2
expect_equivalent(ef$vectors$arrows[1,], c(1,0)) # Dim1 is parallel to Al

## MDSaddpoints puts points to (almost) correct positions
adds <- MDSaddpoints(m, as.matrix(d))
expect_equivalent(adds$points, m$points, tolerance = 1e-3) # 0.00014 in github

## check subsetting of distances needed in MDSaddpoints
d <- vegdist(varespec)
## designdist2 & rectangular subset of d
expect_silent(d2 <- designdist2(varespec[1:10,], varespec[11:24,],
                                 "(A+B-2*J)/(A+B)", "minimum"))
expect_equivalent(dist2xy(d, 11:24), d2)
## subset of symmetric distances
pick <- sample(nrow(varespec), 6, replace = FALSE)
expect_equivalent(dist2xy(d, pick, "xx"), vegdist(varespec[sort(pick),]))



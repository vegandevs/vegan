### tests for running and handling NMDS basically metaMDS -> monoMDS

data(varespec, varechem)
expect_silent(m <- metaMDS(varespec, trymax=1, trace=3))
expect_silent(d <- metaMDSredist(m))
## monoMDS options
init <- wcmdscale(d) # starting position
expect_silent(mo <- monoMDS(d, init, model = "local", weakties = FALSE))
expect_silent(monoMDS(d, init, model = "linear", stress = 2))
expect_silent(monoMDS(d, init, model = "hybrid", scaling = FALSE, pc = FALSE))
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

## check engine and stepacross
data(dune) # varespec has no sites that 'noshare' species
expect_silent(metaMDS(dune, trymax = 0, trace = 0, noshare = TRUE,
                      engine = MASS::isoMDS))

## Rotate to external variable
ef0 <- envfit(m ~ Al + Humdepth + Ca, varechem)
expect_silent(mAl <- MDSrotate(m, varechem$Al))
ef <- envfit(mAl ~ Al + Humdepth + Ca, varechem)
expect_equal(ef$vectors$r, ef0$vectors$r) # unchanged R2
expect_equivalent(ef$vectors$arrows[1,], c(1,0)) # Dim1 is parallel to Al

## MDSaddpoints puts points to (almost) correct positions
adds <- MDSaddpoints(m, as.matrix(d))
expect_equivalent(adds$points, m$points, tolerance = 1e-4) # 4.4e-5 in test

## check subsetting of distances needed in MDSaddpoints
d <- vegdist(varespec)
## designdist2 & rectangular subset of d
expect_silent(d2 <- designdist2(varespec[1:10,], varespec[11:24,],
                                 "(A+B-2*J)/(A+B)", "minimum"))
expect_equivalent(dist2xy(d, 11:24), d2)
## subset of symmetric distances
pick <- sample(nrow(varespec), 6, replace = FALSE)
expect_equivalent(dist2xy(d, pick, "xx"), vegdist(varespec[sort(pick),]))



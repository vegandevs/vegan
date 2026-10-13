### tests for running and handling NMDS basically metaMDS -> monoMDS

data(varespec, varechem)
expect_silent(m <- metaMDS(varespec, trymax=1, trace=3))
expect_silent(d <- metaMDSredist(m))
## monoMDS options
init <- wcmdscale(d) # starting position
expect_silent(monoMDS(d, init, model = "local", weakties = FALSE))
expect_silent(monoMDS(d, init, model = "linear", stress = 2))
expect_silent(monoMDS(d, init, model = "hybrid", scaling = FALSE, pc = FALSE))
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




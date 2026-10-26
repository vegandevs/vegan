### test scores, including sppscores<-, rda scaling, tidy=TRUE

data(dune)

## sppscores<-. mono/metaMDS in test_NMDS

d <- vegdist(dune, "euclidean")
expect_silent(mdb <- pco(d))
expect_silent(mcp <- capscale(dune ~ 1))
expect_silent(mwc <- wcmdscale(d, eig=TRUE))
expect_silent(mx <- wcmdscale(d)) # only coordinates
expect_silent(mrd <- pca(dune))
## capscale vs. rda
expect_silent(pro <- procrustes(mcp, mrd))
expect_equal(abs(pro$rotation), diag(2))
expect_equal(pro$scale, 1)
expect_equal(pro$ss, 0)
expect_silent(sppscores(mcp) <- dune)
expect_silent(pro <- procrustes(mcp, mrd, scores = "species"))
expect_equal(pro$ss, 0)
expect_equal(abs(pro$rotation), diag(2))
expect_equal(pro$scale, 1)
## dbrda vs. rda
expect_silent(pro <- procrustes(mdb, mrd))
expect_equal(abs(pro$rotation), diag(2))
expect_equal(pro$scale, 1)
expect_equal(pro$ss, 0)
expect_silent(sppscores(mdb) <- dune)
expect_silent(pro <- procrustes(mdb, mrd, scores = "species"))
expect_equal(pro$ss, 0)
expect_equal(abs(pro$rotation), diag(2))
expect_equal(pro$scale, 1)
## wcmdscale vs. rda
const <- sqrt(mrd$tot.chi * (nobs(mrd) - 1))
pro <- procrustes(mwc,
                  scores(mrd, display="sites", const=const, scaling="sites"),
                  choices = 1:2)
expect_equal(abs(pro$rotation), diag(2))
expect_equal(pro$scale, 1)
expect_equal(pro$ss, 0)
expect_silent(sppscores(mwc) <- dune)
expect_silent(
    pro <- procrustes(scores(mwc, dis="sp", scaling="sites", choices=1:2),
                   scores(mrd, dis="sp", scaling="sites", const=1)))
expect_equal(pro$ss, 0)
expect_equal(abs(pro$rotation), diag(2))
expect_equal(pro$scale, 1)
## scores.wmdscale vs. scores.default
expect_equal(scores(mwc), scores(mx))
expect_equal(scores(mwc, choices=1:3), scores(mx, choices=1:3))

## tidy scores
data(dune.env)
## rda
mrd <- rda(dune ~ A1 + Moisture, dune.env)
expect_silent(sco <- scores(mrd, display="all"))
expect_equal(dim(tsco <- scores(mrd, display="all", tidy = TRUE)),
             c(sum(sapply(sco, nrow)), 4))
expect_equal(colnames(tsco), c("score", "label", colnames(sco[[1]])))
## cca: tidy scores have 'weight'
mca <- cca(dune ~ A1 + Moisture, dune.env)
expect_silent(sco <- scores(mca, display="all"))
expect_equal(dim(tsco <- scores(mca, display="all", tidy = TRUE)),
             c(sum(sapply(sco, nrow)), 5))
expect_equal(colnames(tsco), c("score", "label", colnames(sco[[1]]), "weight"))
## dbrda
mdb <- dbrda(dune ~ A1 + Moisture, dune.env, distance="bray")
expect_silent(sco <- scores(mdb, display="all"))
expect_equal(dim(tsco <- scores(mdb, display="all", tidy = TRUE)),
             c(sum(sapply(sco, nrow)), 4))
expect_equal(colnames(tsco), c("score", "label", colnames(sco[[1]])))
expect_length(scores(mdb, display = "species"), 0)
expect_silent(sppscores(mdb) <- decostand(dune, "hellinger"))
expect_equal(dim(scores(mdb, display="species")), c(ncol(dune), 2))
## wcmdscale: already have mwc, tidy scores have 'weight'
expect_silent(sco <- scores(mwc, choices=1:2, display=c("sites","species")))
expect_equal(dim(tsco <- scores(mwc, choices = 1:2, tidy = TRUE)),
             c(sum(sapply(sco, nrow)), 5))
expect_equal(colnames(tsco), c("score", "label", colnames(sco[[1]]), "weight"))
## influence statistics for constrained ordination
if (packageVersion("vegan") >= "2.8") {
    expect_silent(infl <- influence(mrd, addFit = TRUE))
    items <- c("hatvalues", "rstandard", "rstudent", "cooks.distance",
               "fitted")
    expect_equal(names(infl), items)
    expect_silent(infl <- influence(mrd, addFit = TRUE, tidy = TRUE))
    expect_equal(dim(infl), c(prod(dim(dune)), 7))
    expect_equal(colnames(infl), c("site", "species", items))
    expect_silent(infl <- influence(mrd, addFit=FALSE, type="canoco",
                                    tidy=TRUE)) # no fitted in type="canoco"
    expect_equal(dim(infl), c(mrd$CCA$rank * nrow(dune), 6))
    expect_equal(colnames(infl), c("site", "axis", items[-5])) # no fitted
}
## scores.default: mx was plain score matrix from wcmdscale(..., eig=FALSE)
expect_equal(dim(scores(mx, choices=1:3, tidy=TRUE)), c(nrow(mx), 5))
## scores.wcmdscale (with eig = TRUE) and added sppscores
expect_silent(tsco <- scores(mwc, choices=1:2, tidy=TRUE))
expect_equal(dim(tsco), c(sum(dim(dune)), 5))
expect_equal(colnames(tsco), c("score","label", colnames(mwc$points[,1:2]),
                               "weight"))
## decorana
expect_silent(mde <- decorana(dune, iweigh = 1)) # downweight: changes weights
expect_silent(sco <- scores(mde, display="both", choices=1:2))
expect_silent(tsco <- scores(mde, choices=1:2, tidy=TRUE))
expect_equal(dim(tsco), c(sum(dim(dune)), 5))
expect_equal(colnames(tsco), c("score","label",colnames(sco[[1]]), "weight"))

## metaMDS & monoMDS (scores.monoMDS delegates to scores.metaMDS)
expect_silent(mds <- metaMDS(dune, trace=0, trymax=0))
expect_silent(sco <- scores(mds, shrink=TRUE))
expect_silent(tsco <- scores(mds, shrink=TRUE, tidy=TRUE))
expect_equal(dim(tsco), c(sum(sapply(sco, nrow)), 4))
expect_equal(colnames(tsco), c("score","label", colnames(sco[[1]])))
## with PCoA init, monoMDS should be equal metaMDS(..., trymax=0)
## except for scaling and possible reflexion of axes
d <- vegdist(dune)
expect_silent(mon <- monoMDS(d, wcmdscale(d)))
expect_silent(sppscores(mon) <- dune)
expect_equal(procrustes(scores(mds), scores(mon))$ss, 0)
expect_equal(procrustes(scores(mon, display="species", shrink=TRUE),
                        scores(mds, display="species", shrink=TRUE))$ss, 0)

## envfit
expect_silent(ef <- envfit(mds ~ ., dune.env))
expect_silent(sco <- scores(ef, arrow.mul = 1))
expect_silent(tsco <- scores(ef, arrow.mul = 1, tidy=TRUE))
expect_equal(dim(tsco), c(sum(sapply(sco, nrow)), 4))
## tidy names are inconsistent with other vegan functions, but
## consistent with ggvegan::tidy
expect_equal(colnames(tsco), c("label","type", colnames(sco[[1]])))

## rarecurve: not a scores method, but can return tidy=TRUE
data(mite)
i <- sample(nrow(mite), 5)
expect_equal(dim(rarecurve(mite[i,], tidy=TRUE)), c(sum(mite[i,]), 3))

## non-ordination scores methods

## betadisper
expect_silent(m <- betadisper(vegdist(dune), dune.env$Management))
expect_silent(sco <- scores(m))

## betadiver
expect_silent(m <- betadiver(mite, method=NA))
Expect_silent(sco <- scores(m, triangular=TRUE))
## all points are within unit-side equilateral triangle, EPS is a
## buffer for numerical accuracy of points at the triangle border.
EPS <- 100 * sqrt(.Machine$double.eps)
triangle <- matrix(c(-EPS, -EPS, 1+EPS, 0-EPS, 0.5, sqrt(0.75)+EPS),
                   nrow=3, ncol=2, byrow=TRUE)
expect_equal(.C(vegan:::pnpoly,
                3L, as.double(triangle[,1]), as.double(triangle[,2]),
                as.integer(nrow(sco)), as.double(sco[,1]), as.double(sco[,2]),
                inpoly = integer(nrow(sco)))$inpoly, # inpoly=1 in triangle
             rep(1, nrow(sco)))

## hclust
cl <- hclust(vegdist(dune))
## internal nodes are at fusion heights
expect_equal(scores(cl, "internal")[,"height"], cl$height)
## x-coordinates of leaves (terminal nodes) are ordered correctly
expect_equal(order(scores(cl, "leaves")[,"x"]), cl$order)

## ordihull: not really intended for normal use, but ordispider uses
## this to draw lines from the group centroid to points on the hull.
mca <- cca(dune)
plot(mca, display="sites")
expect_silent(pl <- ordihull(mca, dune.env$Use, draw="none"))
expect_silent(spi <- ordispider(pl)) # calls scores.ordihull(pl)
expect_equal(rownames(scores(pl)), rownames(spi))

## ordiplot: plot functions return ordiplot object, and its elements
## have scaling options which are retained in scores.ordiplot.
pl <- plot(mca, scaling = "sites", hill = TRUE)
expect_equivalent(scores(pl, display="sites"),
                  scores(mca, display="sites", scaling="sites", hill=TRUE))

## pcnm
data(mite.xy)
expect_equal(dim(scores(pcnm(vegdist(mite.xy, "euclidean")), choices=1:5)),
             c(nrow(mite.xy), 5))








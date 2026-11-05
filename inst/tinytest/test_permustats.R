### permustats cases
data(dune, dune.env)
data(mite, mite.env, mite.xy)
data(sipoo)

## anosim
perm <- shuffleSet(nrow(dune), 49)
m <- anosim(dune, dune.env$Management, permutations = perm)
expect_silent(p <- permustats(m))
expect_equivalent(summary(p)[c("statistic","p")],
                  m[c("statistic","signif")])
## c-method
m1 <- anosim(dune, dune.env$Management, permutations = perm[1:25,])
m2 <- anosim(dune, dune.env$Management, permutations = perm[26:49,])
expect_equal(c(permustats(m1), permustats(m2)),
             permustats(m))

## anova.cca & permutest.cca
perm <- shuffleSet(nrow(dune), 49)
m <- cca(dune ~ A1 + Moisture, dune.env)
ano <- anova(m, permutations = perm)
expect_silent(p <- permustats(ano))
expect_equivalent(unlist(summary(p)["statistic"]), ano[1,"F"])
expect_equivalent(unlist(summary(p)["p"]), ano[1,"Pr(>F)"])
ano <- anova(m, by = "onedf", permutations = perm)
expect_silent(p <- permustats(ano))
expect_equivalent(unlist(summary(p)["statistic"]),
                  ano[seq_len(m$CCA$rank),"F"])
expect_equivalent(unlist(summary(p)["p"]),
                  ano[seq_len(m$CCA$rank),"Pr(>F)"])
## density, qqnorm & boxplot: numerical tests
expect_error(density(p), "more than one statistic")
expect_silent(density(p, which=2))
expect_error(qqnorm(p), "more than one statistic")
expect_silent(qq <- qqnorm(p, observed=FALSE, which=2))
expect_equal(qq$y, sort(p$permutations[,2]))
expect_silent(bp <- boxplot(p, range = 0))
## boxplot centres so that statistic==0
expect_equal(sweep(bp$stats, 2, p$statistic, "+"),
             apply(p$permutations, 2, fivenum))
## c-method
an1 <- anova(m, by="onedf", permutations = perm[1:24,])
an2 <- anova(m, by="onedf", permutations = perm[25:49,])
expect_equal(c(permustats(an1), permustats(an2)),
             permustats(ano))
## permutest.cca
ano <- permutest(m, by="term")
expect_silent(p <- permustats(ano))
## P-value is only evaluated within print.permustest.cca
expect_equivalent(unlist(summary(p)$statistic), ano$F.0)

## oecosimu family
nqsw <- nullmodel(sipoo, "quasiswap")
mods <- simulate(nqsw, 49)
expect_silent(oe <- oecosimu(mods, nestedn0))
expect_silent(p <- permustats(oe))
expect_equivalent(summary(p)[c("statistic", "p")],
                  oe$oecosimu[c("statistic", "pval")])
## subsetting simmat fails mods[,,1:10] returns 3D array, but drops
## class and other simmat attributes - therefore no c() tests here.

## permutest.betadisper
d <- vegdist(dune)
m <- betadisper(d, dune.env$Management, bias.adjust = TRUE)
expect_silent(p <- permutest(m, permutations = 99))
expect_silent(sum <- summary(permustats(p)))
expect_equivalent(unlist(sum[c("statistic", "p")]),
                  p$tab[1,c("F","Pr(>F)")])

## protest
mod1 <- cca(dune ~ Moisture + Management, dune.env)
mod2 <- cca(dune ~ A1 + Moisture + Management + Use, dune.env)
m <- protest(mod1, mod2, permutations = 99)
expect_silent(sum <- summary(permustats(m)))
expect_equivalent(sum[c("statistic", "p")], m[c("t0","signif")])

## the following methods use summary.permustats as summary

## mantel & mantel.partial
dis <- vegdist(sqrt(mite), "hellinger")
gdis <- vegdist(mite.xy, "euclid")
edis <- vegdist(mite.env[,1:2], "euclid")
m <- mantel(dis, edis, permutations = 99)
expect_silent(sum <- summary(m))
expect_equivalent(sum[c("statistic","p")], m[c("statistic","signif")])
m <- mantel.partial(dis, edis, gdis, permutations = 99)
expect_silent(sum <- summary(m))
expect_equivalent(sum[c("statistic","p")], m[c("statistic","signif")])

## mrpp
m <- mrpp(dune, dune.env$Management, weight.type = 3, permutations = 99)
expect_silent(p <- summary(m))
expect_equivalent(p[c("statistic", "p")], m[c("delta","Pvalue")])

## ordiareatest needs plot on which area is evaluated
mod <- cca(dune)
plot(mod, display="sites")
m <- ordiareatest(mod, dune.env$Management, kind = "ellipse", permutations = 99)
expect_silent(sum <- summary(m))
expect_equivalent(sum[c("statistic", "p")], m[c("areas","pvalues")])

## the following methods return no permutation data

if (packageVersion("vegan") >= "2.8") {
## parameteric anova.cca
mod <- rda(dune ~ Management + Moisture, dune.env)
ano <- anova(mod, by="margin", test = "parametric")
expect_error(permustats(ano), "no permutation data available")
}

## CCorA
m <- CCorA(mite, mite.env, permutations = 9)
expect_error(permustats(m), "no permutation data available")

## envfit, vectorfit, factorfit
ef <- envfit(mod ~ Moisture + A1, dune.env)
expect_error(permustats(ef), "no permutation data available")
expect_error(permustats(ef$vectors), "no permutation data available")
expect_error(permustats(ef$factors), "no permutation data available")

## mso
mod <- cca(mite ~ Topo + Shrub + WatrCont, mite.env)
expect_silent(m <- mso(mod, mite.xy, permutations = 19))
expect_error(permustats(m), "no permutation data available")

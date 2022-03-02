# Test data
# data(varespec)
testdata <- matrix(round(runif(1000, 0, 100)), nrow=20)
testdata <- testdata - 50
testdata[testdata < 0] <- 0
rownames(testdata) <- paste0("row", seq_len(nrow(testdata)))
colnames(testdata) <- paste0("col", seq_len(ncol(testdata)))

# Calculates relative abundance table
relative <- vegan::decostand(testdata, "total")

# Count and relative data with pseudocount
testdata.with.pseudo <- testdata + 1
relative.with.pseudo <- vegan::decostand(testdata+1, "total")

# aitchison equals to CLR + Euclid (pseudocount is necessary with clr)
a1 <- vegan::vegdist(testdata+1, method = "aitchison")
a2 <- vegan::vegdist(vegan::decostand(testdata+1, "clr"), method = "euclidean")
max(abs(a1-a2)) < 1e-6 # Tolerance

# Robust aitchison equals to rCLR + Euclid
# and works without pseudocount
a1 <- vegan::vegdist(testdata, method = "robust.aitchison")
a2 <- vegan::vegdist(vegan::decostand(testdata, "rclr"), method = "euclidean")
max(abs(a1-a2)) < 1e-6 # Tolerance

# Robust aitchison and aitchison are equal when there are no zeroes
a1 <- vegan::vegdist(testdata.with.pseudo, method = "robust.aitchison")
a2 <- vegan::vegdist(testdata.with.pseudo, method = "aitchison")
max(abs(a1-a2)) < 1e-6 # Tolerance

# It is possible to pass pseudocount as a function argument to vegan::decostand
a1 <- vegan::vegdist(testdata, method = "aitchison", pseudocount=1)
a2 <- vegan::vegdist(testdata+1, method = "aitchison")
max(abs(a1-a2)) < 1e-6 # Tolerance


# Compare the outcomes with an external package that also provides compositional transformations
# Adding these would demand adding Suggested packages in DESCRIPTION; skipped for now but can be
# useful for manual testing.
#skip <- TRUE
#if (!skip) {
#
#    sum(compositions::clr(testdata.with.pseudo) - vegan::decostand(testdata.with.pseudo, "clr")) < 1e-6
#    sum(compositions::clr(testdata.with.pseudo) - vegan::decostand(testdata.with.pseudo, "clr", MARGIN=1)) < 1e-6
#    sum(t(compositions::clr(t(testdata.with.pseudo))) - vegan::decostand(testdata.with.pseudo, "clr", MARGIN=2)) < 1e-6
#    sum(rgr::clr(testdata.with.pseudo) - vegan::decostand(testdata.with.pseudo, "clr"))<1e-6#
#
#    sum(compositions::alr(testdata.with.pseudo, ivar=1) - vegan::decostand(testdata.with.pseudo, "alr")) < 1e-6
#    sum(compositions::alr(testdata.with.pseudo, ivar=1) - vegan::decostand(testdata.with.pseudo, "alr", MARGIN=1)) < 1e-6
#    sum(t(compositions::alr(t(testdata.with.pseudo), ivar=1)) - vegan::decostand(testdata.with.pseudo, "alr", MARGIN=2)) < 1e-6
#    sum(rgr::alr(testdata.with.pseudo, j=1) - vegan::decostand(testdata.with.pseudo, "alr", reference=1))<1e-6#
#
#}

# --------------------------------------------------------------------------------------------------------------


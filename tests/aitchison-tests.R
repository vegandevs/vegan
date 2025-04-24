###<--- BEGIN TESTS --->
suppressPackageStartupMessages(require(vegan))

# Test data
# data(varespec)
set.seed(252)
testdata <- matrix(round(runif(1000, 0, 100)), nrow = 20)
testdata <- testdata - 50
testdata[testdata < 0] <- 0
rownames(testdata) <- paste0("row", seq_len(nrow(testdata)))
colnames(testdata) <- paste0("col", seq_len(ncol(testdata)))

# Calculates relative abundance table
relative <- decostand(testdata, "total")

# Count and relative data with pseudocount
testdata.with.pseudo <- testdata + 1
relative.with.pseudo <- decostand(testdata + 1, "total")

# NAs correctly assigned
a1 <- testdata
a2 <- decostand(a1, "rclr", na.rm = FALSE)
a3 <- decostand(a1, "rclr", na.rm = TRUE)
all(is.na(a2[a1 == 0]))
all(!is.na(a3[a1 == 0]))

# aitchison equals to CLR + Euclid (pseudocount is necessary with clr)
a1 <- vegdist(testdata + 1, method = "aitchison", na.rm = TRUE)
a2 <- vegdist(decostand(testdata + 1, "clr"), method = "euclidean", na.rm = TRUE)
max(abs(a1-a2)) < 1e-6 # Tolerance

# Robust aitchison equals to rCLR + Euclid when na.rm = TRUE (matrix completion applied)
# and works without pseudocount
a1 <- vegdist(testdata, method = "robust.aitchison", na.rm = TRUE)
a2 <- vegdist(decostand(testdata, "rclr"), method = "euclidean", na.rm = TRUE)
max(abs(a1-a2)) < 1e-6 # Tolerance

# Robust aitchison and aitchison are equal when there are no zeroes
a1 <- vegdist(testdata.with.pseudo, method = "robust.aitchison")
a2 <- vegdist(testdata.with.pseudo, method = "aitchison")
max(abs(a1 - a2)) < 1e-6 # Tolerance

# It is possible to pass pseudocount as a function argument to decostand
a1 <- vegdist(testdata, method = "aitchison", pseudocount = 1)
a2 <- vegdist(testdata + 1, method = "aitchison")
max(abs(a1 - a2)) < 1e-6 # Tolerance

# Robust rclr+Euclid should give same result than robust.aitchison
d1 <- vegdist(decostand(testdata, "rclr"), "euclidean", na.rm = TRUE)
d2 <- vegdist(testdata, "robust.aitchison", na.rm = TRUE)
max(abs(d1 - d2)) < 1e-6 # Tolerance

# Robust rclr+Euclid should give same result than robust.aitchison
d1 <- vegdist(decostand(testdata + 1, "rclr"), "euclidean", na.rm = FALSE)
d2 <- vegdist(testdata + 1, "robust.aitchison", na.rm = FALSE)
max(abs(d1 - d2)) < 1e-6 # Tolerance

# Compare the outcomes with an external package that also provides compositional transformations
# Adding these would demand adding Suggested packages in DESCRIPTION
# skipped for now but can be useful for manual testing.
#skip <- TRUE
#if (!skip) {
#
#    sum(compositions::clr(testdata.with.pseudo) - decostand(testdata.with.pseudo, "clr")) < 1e-6
#    sum(compositions::clr(testdata.with.pseudo) - decostand(testdata.with.pseudo, "clr", MARGIN = 1)) < 1e-6
#    sum(t(compositions::clr(t(testdata.with.pseudo))) - decostand(testdata.with.pseudo, "clr", MARGIN = 2)) < 1e-6
#    sum(rgr::clr(testdata.with.pseudo) - decostand(testdata.with.pseudo, "clr"))<1e-6#
#
#    sum(compositions::alr(testdata.with.pseudo, ivar = 1) - decostand(testdata.with.pseudo, "alr")) < 1e-6
#    sum(compositions::alr(testdata.with.pseudo, ivar = 1) - decostand(testdata.with.pseudo, "alr", MARGIN = 1)) < 1e-6
#    sum(t(compositions::alr(t(testdata.with.pseudo), ivar = 1)) - decostand(testdata.with.pseudo, "alr", MARGIN = 2)) < 1e-6
#    sum(rgr::alr(testdata.with.pseudo, j = 1) - decostand(testdata.with.pseudo, "alr", reference = 1))<1e-6#
#
#}

# --------------------------------------------------------------------------------------------------------------

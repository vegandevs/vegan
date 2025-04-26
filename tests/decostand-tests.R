###<--- BEGIN TESTS --->
suppressPackageStartupMessages(require(vegan))

############################### CLR ####################################
# Calculates clr-transformation. Should be equal.

# Test data
data(varespec)
set.seed(42)
testdata <- matrix(round(runif(1000, 0, 100)), nrow = 20)
testdata <- testdata - 50
testdata[testdata < 0] <- 0
rownames(testdata) <- paste0("row", seq_len(nrow(testdata)))
colnames(testdata) <- paste0("col", seq_len(ncol(testdata)))
testdata.with.pseudo <- testdata + 1

# Calculates relative abundance table
relative <- decostand(testdata, "total")
relative.with.pseudo <- decostand(testdata + 1, "total")

# CLR data
x.clr  <- decostand(testdata + 1, method = "clr")
x.rclr <- decostand(testdata, method = "rclr")
x.clr.pseudo <- decostand(testdata, method = "clr", pseudocount = 1)

max(abs(x.clr - x.clr.pseudo)) < 1e-6
max(abs(decostand(testdata + 1, method = "clr", pseudocount = 0) - x.clr.pseudo)) < 1e-6

# Tests clr
alt.clr <- t(apply(as.matrix(relative.with.pseudo), 1, FUN = function(x){log(x) - mean(log(x))}))
max(abs(x.clr - alt.clr)) < 1e-6

# Test that NAs are handled as expected in CLR
x <- testdata
set.seed(24)
x[sample(prod(dim(x)), 50)] <- NA # insert some NAs in the data
# NAs in the original data remain NAs in the returned data
all(is.na(decostand(x, "clr", na.rm = FALSE, pseudocount = 1)[is.na(x)])) == TRUE # NAs
# For the other (non-NA) values, we get non-NA values back
any(is.na(decostand(x, "clr", na.rm = FALSE, pseudocount = 1)[!is.na(x)])) == FALSE
any(is.na(decostand(x, "clr", MARGIN = 2, na.rm = FALSE, pseudocount = 1)[!is.na(x)])) == FALSE 
# Results match for the non-NA values always (with tolerance 1e-6)
inds <- !is.na(x) # Non-NA values
max(abs(decostand(x, "clr", na.rm = FALSE, pseudocount = 1)[inds] -
        decostand(x, "clr", na.rm = TRUE, pseudocount = 1)[inds])) < 1e-6
# For the other (non-NA) values, we get non-NA values back
any(is.na(decostand(x, "alr", na.rm = FALSE, pseudocount = 1)[!is.na(x[, -1] - x[, 1])])) == FALSE
any(is.na(decostand(x, "alr", na.rm = FALSE, pseudocount = 1, reference=3)[!is.na(x[, -3] - x[, 3])])) == FALSE
# Works correctly also with other MARGIN
inds <- !is.na(x) # Non-NA values
max(abs(decostand(x, "clr", MARGIN = 2, na.rm = FALSE, pseudocount = 1)[inds] -
        decostand(x, "clr", MARGIN = 2, na.rm = TRUE, pseudocount = 1)[inds])) < 1e-6
# Results match for the non-NA values always (with tolerance 1e-6)
inds <- !is.na(x) # Non-NA values
max(abs(na.omit(decostand(x, "alr", na.rm = FALSE, pseudocount = 1)[inds] -
                decostand(x, "alr", na.rm = TRUE, pseudocount = 1)[inds]))) < 1e-6

# Test that NAs are handled as expected in ALR
set.seed(4354353)
x <- testdata
x[sample(prod(dim(x)), 50)] <- NA
x[4, c(2, 10)] <- NA  # insert some NAs in the data
# NAs in the output share NAs with the reference vector
all(is.na(decostand(x, "alr", na.rm = FALSE, pseudocount = 1, reference = 4))[which(is.na(x[,4])),])
# Output vector has same NAs than the original vector and reference vector
all(is.na(decostand(x, "alr", na.rm = FALSE, pseudocount = 1, reference = 4)[,2]) == (is.na(x[,4]) | is.na(x[,2])))
# No NAs after removing them
!any(is.na(decostand(x, "alr", na.rm = TRUE, pseudocount = 1, reference = 4)))
# All NAs are replaced by zero
all((rowMeans(decostand(x, "alr", na.rm = TRUE, pseudocount = 1, reference = 4) == 0) == 1) == is.na(x[,4]))

# Expect that error does not occur
decostand(testdata, method = "rclr")
decostand(testdata, method = "clr", pseudocount = 1)
# Expect error
#class(try(decostand(testdata, method = "clr"))) == "try-error"
#class(try(decostand(testdata, method = "clr", pseudocount = 0))) == "try-error"

# Tests that clr robust gives values that are approximately same if only 
# one value per sample are changed to zero
# Adds pseudocount
test <- testdata + 1
test2 <- test
test2[,1] <- 0

# clr robust transformations
test <- decostand(test, method = "rclr")
test2 <- decostand(test2, method = "rclr")

# Removes first cols
test <- test[, -1]
test2 <- test2[, -1]

# Expect high correlation
cor(unlist(test), unlist(test2)) > 0.99

# Compare rclr with imputation to without imputation + matrix completion
# rclr transformation with matrix completion for the 0/NA entries
x1 <- decostand(varespec, method = "rclr", impute = TRUE)

# rclr transformation with no matrix completion for the 0/NA entries
x2 <- decostand(varespec, method = "rclr", impute = FALSE)

# Matrix completion
x2c <- optspace(x2, ropt = 3, niter = 5, tol = 1e-5, verbose = FALSE)$M
all(as.matrix(x1) == as.matrix(x2c))

x2c <- optspace(x2, niter = 5, tol = 1e-5, verbose = FALSE)$M
all(as.matrix(x1) == as.matrix(x2c))

############################# NAMES ####################################

# Tests that dimensins have correct names
all(rownames(decostand(testdata + 1, method = "clr")) == rownames(testdata))
all(colnames(decostand(testdata, method = "clr", pseudocount = 1)) == colnames(testdata))
all(rownames(decostand(testdata, method = "rclr")) == rownames(testdata))
all(colnames(decostand(testdata, method = "rclr")) == colnames(testdata))

########################################################################

# Count vs. Relative data

# CLR is identical with count and relative data
a1 <- decostand(testdata.with.pseudo, method = "clr")
a2 <- decostand(relative.with.pseudo, method = "clr")
max(abs(a1 - a2)) < 1e-6 # Tolerance

# rCLR is identical with count and relative data
a1 <- decostand(testdata.with.pseudo, method = "rclr")
a2 <- decostand(relative.with.pseudo, method = "rclr")
max(abs(a1 - a2)) < 1e-6 # Tolerance

# ALR is identical with count and relative data
a1 <- decostand(testdata.with.pseudo, method = "alr")
a2 <- decostand(relative.with.pseudo, method = "alr")
max(abs(a1 - a2)) < 1e-6 # Tolerance

####### # ALR transformation drops one feature ################
ncol(decostand(testdata.with.pseudo, "alr")) == ncol(testdata.with.pseudo) - 1

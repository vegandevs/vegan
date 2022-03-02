        
############################### CLR ####################################
# Calculates clr-transformation. Should be equal.

# Test data
# data(varespec)
testdata <- matrix(round(runif(1000, 0, 100)), nrow=20)
testdata <- testdata - 50
testdata[testdata < 0] <- 0
rownames(testdata) <- paste0("row", seq_len(nrow(testdata)))
colnames(testdata) <- paste0("col", seq_len(ncol(testdata)))
testdata.with.pseudo <- testdata + 1

# Calculates relative abundance table
relative <- vegan::decostand(testdata, "total")
relative.with.pseudo <- vegan::decostand(testdata+1, "total")

# CLR data
x.clr <- vegan::decostand(testdata+1, method = "clr")
x.rclr <- vegan::decostand(testdata, method = "rclr")
x.clr.pseudo <- vegan::decostand(testdata, method = "clr", pseudocount=1)

max(abs(x.clr - x.clr.pseudo))<1e-6
max(abs(vegan::decostand(testdata+1, method = "clr", pseudocount=0)-x.clr.pseudo))<1e-6

# Tests clr
alt.clr <- t(apply(as.matrix(relative.with.pseudo), 1, FUN=function(x){
log(x) - mean(log(x))}))
max(abs(x.clr-alt.clr)) < 1e-6      
all((x.rclr==0) == (testdata==0))

# Expect that error does not occur
vegan::decostand(testdata, method = "rclr")
vegan::decostand(testdata, method = "clr", pseudocount=1)
# Expect error
class(try(vegan::decostand(testdata, method = "clr")))=="try-error"
class(try(vegan::decostand(testdata, method = "clr", pseudocount=0)))=="try-error"

# Tests that clr robust gives values that are approximately same if only 
# one value per sample are changed to zero
# Adds pseudocount
test <- testdata+1 
test2 <- test; test2[,1] <- 0

# clr robust transformations
test <- vegan::decostand(test, method = "rclr")
test2 <- vegan::decostand(test2, method = "rclr")

# Removes first cols
test <- test[, -1]
test2 <- test2[, -1]

# Expect high correlation
cor(unlist(test), unlist(test2)) > 0.99

############################# NAMES ####################################

# Tests that dimensins have correct names
all(rownames(vegan::decostand(testdata+1, method = "clr")) == rownames(testdata))
all(colnames(vegan::decostand(testdata, method = "clr", pseudocount=1)) == colnames(testdata))
all(rownames(vegan::decostand(testdata, method = "rclr")) == rownames(testdata))
all(colnames(vegan::decostand(testdata, method = "rclr")) == colnames(testdata))

########################################################################

# Count vs. Relative data

# CLR is identical with count and relative data
a1 <- vegan::decostand(testdata.with.pseudo, method = "clr")
a2 <- vegan::decostand(relative.with.pseudo, method = "clr")
max(abs(a1-a2)) < 1e-6 # Tolerance

# rCLR is identical with count and relative data
a1 <- vegan::decostand(testdata.with.pseudo, method = "rclr")
a2 <- vegan::decostand(relative.with.pseudo, method = "rclr")
max(abs(a1-a2)) < 1e-6 # Tolerance

# ALR is identical with count and relative data
a1 <- vegan::decostand(testdata.with.pseudo, method = "alr")
a2 <- vegan::decostand(relative.with.pseudo, method = "alr")
max(abs(a1-a2)) < 1e-6 # Tolerance

####### # ALR transformation drops one feature ################
ncol(vegan::decostand(testdata.with.pseudo, "alr")) == ncol(testdata.with.pseudo)-1







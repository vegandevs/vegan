### oecosimu-tests: unit tests for vegan functions

### This file contains unit tests for oecosimu, permatfull and related
### functions. This file is run in R CMD check and its results are
### compared against previously saved results in
### oecosimu-tests.Rout.save. If you change tests, you must generate
### new oecosimu-tests.Rout.save in this directory.

### The current plan is that tests/ are not included in the CRAN
### release, but only in the development version of vegan in R-Forge.

### The tests here are not intended for human reading. The tests need
### not be ecological or biologically meaningful, but they are only
### intended for testing strange arguments, protect against
### regressions and test correctness of results.

### The tests are in a single file and you should clean work space
### after the unit test. You should set random number seed (if needed)
### for comparison against vegan-tests.Rout.save, and you should
### remove the seed after the unit test. If possible, avoid very long
### lasting tests.

### <-- oecosimu/permat specifics -->

### The necessary condition is that fill, row and column sums are honoured. Secondarily we want to keep the simulation results identical when chaning the underlying routines in oecosimu (except experimental method = "backtrack"). The permatfull/swap code can be regarded as more experimental, and we can perhaps there change the simulation results, but not lightly. 

###<--- BEGIN TESTS --->
suppressPackageStartupMessages(require(vegan))
set.seed(4711)

### Prior to vegan 2.2-0 these were commsimulator() tests for binary
### null model, but now commsimulator() is deprecated and calls were
### replaced with corresponding simulate(nullmodel())

## expect fill, rowSums, colSums
expect <- data.frame("r00" = c(TRUE, FALSE, FALSE),
                     "r0" = c(TRUE, TRUE, FALSE),
                     "r1" = c(TRUE, TRUE, FALSE),
                     "r2" = c(TRUE, TRUE, FALSE),
                     "c0" = c(TRUE, FALSE, TRUE),
                     "swap" = c(TRUE, TRUE, TRUE),
                     "tswap" = c(TRUE, TRUE, TRUE),
                     "quasiswap" = c(TRUE, TRUE, TRUE))
margintest <- function(x, fill, rs, cs) {
    c(sum(x) == fill,
      all(rowSums(x) == rs),
      all(colSums(x) == cs))
}

## data must be binary! 'sipoo' is (and not transformed) 
data(sipoo)
fill <- sum(sipoo)
rs <- rowSums(sipoo)
cs <- colSums(sipoo)

for(method in names(expect)) {
    m <- simulate(nullmodel(sipoo, method = method), thin = 100)
    cat("--> method:", method, "\n")
    cat("** margintest:")
    print(all(margintest(m, fill, rs, cs) == expect[, method]))
    vegemite(drop(m)) 
}
## clean
rm(list = ls())

### end commsimulator

### set up permatfull/swap tests

data(mite)
data(mite.env)
x <- as.matrix(subset(mite, mite.env$Topo == "Blanket"))
x <- x[, colSums(x) > 0]
gsum <- sum(x)
rsum <- rowSums(x)
csum <- colSums(x)
fill <- sum(x>0)
rfrq <- rowSums(x>0)
cfrq <- colSums(x>0)

margintest <-
    function(x, gsum, rsum, csum, fill, rfrq, cfrq)
{
    cat("grand sum: ")
    print(all.equal(sum(x), gsum, check.attributes = FALSE))
    cat("row sums:  ")
    print(all.equal(rowSums(x), rsum, check.attributes = FALSE))
    cat("col sums:  ")
    print(all.equal(colSums(x), csum, check.attributes = FALSE))
    cat("fill:      ")
    print(all.equal(sum(x>0), fill, check.attributes = FALSE))
    cat("row freqs: ")
    print(all.equal(rowSums(x>0), rfrq, check.attributes = FALSE))
    cat("col freqs: ")
    print(all.equal(colSums(x>0), cfrq, check.attributes = FALSE))
}

### permatfull

set.seed(4711)
margin <- c("none", "rows", "columns", "both")
shuffle <- c("ind", "samp", "both")

for(mar in margin) {
    for(what in shuffle) {
        cat("\n--> margin", mar, " shuffle", what, "<--\n")
        set.seed(4711)
        m <- permatfull(x, fixedmar = mar, shuffle = what, mtype = "count", 
            times=1)$perm[[1]]
        margintest(m, gsum, rsum, csum, fill, rfrq, cfrq)
        print(m[,1:12])
    }
}

### permatswap
set.seed(4711)
methods <- c("swap", "quasiswap", "swsh", "abuswap")
margins <- c("rows", "columns", "both")
shuffle <- c("samp", "both")

incompatible <- function(method, margin)
{
    (method == "swap" && margin != "both") ||
    (method == "abuswap" && margin == "both") ||
    (method == "quasiswap" && margin != "both") ||
    (method == "swsh" && margin == "both")
}

for(method in methods) {
    for(margin in margins) {
        for(what in shuffle) {
            if (incompatible(method = method, margin = margin))
                next
            cat("\n*** ", method, " ***\n")
            cat("--> margin", margin, " shuffle", what, "<--\n")
            set.seed(4711)
            m <- permatswap(x, method = method, fixedmar = margin, shuffle = what,
                             mtype = "count", thin=100, times=1)$perm[[1]]
            margintest(m, gsum, rsum, csum,fill, rfrq, cfrq)
            print(m[,1:12])
        }
    }
}
### end permatswap

nm=nullmodel(m,"abuswap_c")
sm<-simulate(nm,nsim=100,thin=100,burnin=100)

### clean
rm(list = ls())

## end permatfull/swap

## The following vegan functions depend on *oecosimu*: adipart
## hiersimu multipart raupcrick. The following functions directly
## depend on *commsimulator*: permatfull1 permatswap1.  All these have
## derived and/or method functions. These should not be broken.

## Do not break raupcrick:
set.seed(4711)
data(sipoo)
as.numeric(raupcrick(sipoo, nsimul = 99))
rm(list = ls())
## end raupcrick



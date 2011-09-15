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

### start commsimulator

## expect fill, rowSums, colSums
expect <- data.frame("r00" = c(TRUE, FALSE, FALSE),
                     "r0" = c(TRUE, TRUE, FALSE),
                     "r1" = c(TRUE, TRUE, FALSE),
                     "r2" = c(TRUE, TRUE, FALSE),
                     "c0" = c(TRUE, FALSE, TRUE),
                     "swap" = c(TRUE, TRUE, TRUE),
                     "tswap" = c(TRUE, TRUE, TRUE),
                     "quasi" = c(TRUE, TRUE, TRUE))
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
    m <- commsimulator(sipoo, method = method, thin = 100)
    cat("--> method:", method, "\n")
    cat("** margintest:")
    print(all(margintest(m, fill, rs, cs) == expect[, method]))
    vegemite(m)
}
## clean
rm(list = ls())

### end commsimulator

### oecosimu-tests: unit tests for vegan functions

### This file contains basic unit tests for simulating null
### models. Currently we just test the marginal properties of null
### models using example(commsim).

### We had more extensive tests that also displayed heads of simulated
### matrices (not only the marginal sums), but these were unstable
### when R was compiled as ./configure --disable-long-double because
### some tests used R functions stats::rmultinom() which used long
### doubles and gave different results when long double was not
### available.

### We also had specific permatfull/permatswap tests, but these only
### tested the simple cases and were nothing but an alternative
### interface to commsim nullmodels tested here.

### <-- oecosimu/permat specifics -->

###<--- BEGIN TESTS --->
suppressPackageStartupMessages(require(vegan))
set.seed(4711)
example(commsim)

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

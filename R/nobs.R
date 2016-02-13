### R 2.13.0 introduces nobs() method to get the number of
### observations. This file provides methods for vegan classes.

`nobs.adonis` <- function(object, ...) NROW(object$coef.sites)

`nobs.betadisper` <- function(object, ...) length(object$distances)

`nobs.cca` <- function(object, ...) max(NROW(object$pCCA$u),
                                        NROW(object$CCA$u),
                                        NROW(object$CA$u))

`nobs.CCorA` <- function(object, ...) NROW(object$Cy)

`nobs.decorana` <- function(object, ...) NROW(object$rproj)

`nobs.isomap` <- function(object, ...) NROW(object$points)

`nobs.metaMDS` <- function(object, ...) NROW(object$points)

`nobs.pcnm` <- function(object, ...) NROW(object$vectors)

`nobs.procrustes` <- function(object, ...) NROW(object$X)

`nobs.rad` <- function(object, ...) length(object$y)

`nobs.varpart` <- function(object, ...) object$part$n

`nobs.wcmdscale` <- function(object, ...) NROW(object$points)

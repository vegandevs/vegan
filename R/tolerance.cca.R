##' Species tolerances and sample heterogeneities
##'
##' Function to compute species tolerances and site heterogeneity measures
##' from unimodal ordinations (CCA & CA). Implements Eq 6.47 and 6.48 from
##' the Canoco 4.5 Reference Manual (pages 178-179).
##'
##' @param x object of class \code{"cca"}.
##' @param choices numeric; which ordination axes to compute
##' tolerances  and heterogeneities for. Defaults to axes 1 and 2.
##' @param which character; one of \code{"species"} or \code{"sites"},
##' indicating whether species tolerances or sample heterogeneities
##' respectively are computed.
##' @param scaling numeric or character; the ordination scaling to use.
##' @param useN2 logical; should the bias in the tolerances /
##' heterogeneities be reduced via scaling by Hill's N2?
##' @param ... arguments passed to other methods
##' @return matrix of tolerances/heterogeneities with some additional
##'   attributes: \code{which}, \code{scaling}, and \code{N2}, the latter of which will be \code{NA} if \code{useN2 = FALSE}.
##' @author Gavin L. Simpson
##' @examples
##' data(dune)
##' data(dune.env)
##' mod <- cca(dune ~ ., data = dune.env)
##' tolerance.cca(mod)
##'
tolerance.cca <- function(x, choices = 1:2,
                          which = c("species","sites"),
                          scaling = "species", useN2 = TRUE,
                          hill = FALSE, ...) {
    if(inherits(x, "rda")) {
        stop("tolerances only available for unimodal ordinations")
    }
    if(missing(which)) {
        which <- "species"
    }
    ## zapping epsilon to make approximate 1's into 1's
    ZAP <- sqrt(.Machine$double.eps)
    ## reconstruct species/response matrix Y - up to machine precision!
    Y <- (ordiYbar(x, "initial") * sqrt(x$rowsum %o% x$colsum) +
          x$rowsum %o% x$colsum) * x$grand.total
    which <- match.arg(which)
    siteScrTypes <- if (is.null(x$CCA)) {
                        "sites"
                    } else {
                        "lc"
                    }
    ## Sort out scaling; only for (C)CA so no correlation arg
    scaling <- scalingType(scaling, hill = hill)
    scrs <- scores(x, display = c(siteScrTypes, "species"),
                   choices = choices, scaling = scaling, ...)
    ## compute N2 if useN2 == TRUE & only if
    doN2 <- isTRUE(useN2) && ((which == "species" && abs(scaling) == 2) ||
                              (which == "sites" && abs(scaling) == 1))
    siteScrs <- which(names(scrs) %in% c("sites","constraints"))
    if(isTRUE(all.equal(which, "sites"))) {
        res <- matrix(ncol = length(choices), nrow = nrow(scrs[[siteScrs]]))
        Ytot <- rowSums(Y)
        for (i in seq_len(NROW(res))) {
            XiUk <- apply(scrs[["species"]], 1L, `-`, scrs[[siteScrs]][i,])
            YXiUk <- sweep(XiUk^2, 2L, Y[i,], "*")
            if(any(neg <- YXiUk < 0)) {
                YXiUk[neg] <- 0
            }
            res[i, ] <- sqrt(rowSums(YXiUk) / Ytot[i])
        }
        rownames(res) <- rownames(scrs[[siteScrs]])
        colnames(res) <- colnames(scrs[[siteScrs]])

        if(doN2) {
            y <- sweep(Y, 1, Ytot, "/")^2
            N2 <- 1 / rowSums(y, na.rm = TRUE) ## 1/H
            ## avoid almost-1 for sites with only one spp
            N2[abs(N2-1) < ZAP] <- 1
            ## avoid "negative zeros" form 1 - 1/N2 when N2 ~ 1
            res <- sweep(res, 1, sqrt(pmax(1 - 1/N2, 0)), "/")
        }
    } else {
        res <- matrix(ncol = length(choices), nrow = ncol(Y))
        Ytot <- colSums(Y)
        for (i in seq_len(NROW(res))) {
            XiUk <- apply(scrs[[siteScrs]], 1L, `-`, scrs[["species"]][i,])
            YXiUk <- sweep(XiUk^2, 2L, Y[,i], "*")
            if (any(neg <- YXiUk < 0)) {
                YXiUk[neg] <- 0
            }
            res[i, ] <- sqrt(rowSums(YXiUk) / Ytot[i])
        }
        rownames(res) <- colnames(Y)
        colnames(res) <- colnames(scrs[["species"]])

        if(doN2) {
            y <- sweep(Y, 2, Ytot, "/")^2
            N2 <- 1 / colSums(y, na.rm = TRUE) # 1/H
            ## avoid almost-1 for species present only once
            N2[abs(N2-1) < ZAP] <- 1
            ## avoid "negative zeros" form 1 - 1/N2 when N2 ~ 1
            res <- sweep(res, 1, sqrt(pmax(1 - 1/N2, 0)), "/")
        }
    }
    res[!is.finite(res)] <- 0 # some values can be Inf or NaN but are really 0
    res[res < sqrt(.Machine$double.eps)] <- 0 # almost-zero tolerances should be zero
    class(res) <- c("tolerance.cca", "tolerance","matrix")
    attr(res, "which") <- which
    attr(res, "scaling") <- scaling
    attr(res, "N2") <- NA
    if(doN2) {
        attr(res, "N2") <- N2
    }
    res                                 # return
}

`print.tolerance.cca` <- function(x, ...) {
    cat("\n")
    msg <- ifelse(attr(x, "which") == "species", "Species Tolerance",
                  "Sample Heterogeneity")
    writeLines(msg, sep = "\n\n")
    msg <- paste("Scaling:", attr(x, "scaling"))
    writeLines(strwrap(msg), sep = "\n\n")
    attr(x, "model") <- attr(x, "scaling") <- attr(x, "which") <- attr(x, "N2") <- NULL
    print(unclass(x), ...)
    cat("\n")
}

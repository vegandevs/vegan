##' Species tolerances and sample heterogeneities
##'
##' Function to compute species tolerances and site heterogeneity measures
##' from unimodal ordinations (CCA & CA). Implements Eq 6.47 and 6.48 from
##' the Canoco 4.5 Reference Manual (pages 178-179).
##'
##' @param x object of class \code{"cca"}.
##' @param choices numeric; which ordination axes to compute
##' tolerances   and heterogeneities for. Defaults to axes 1 and 2.
##' @param which character; one of \code{"species"} or \code{"sites"},
##' indicating whether species tolerances or sample heterogeneities
##' respectively are computed.
##' @param scaling numeric; the ordination scaling to use.
##' @param useN2 logical; should the bias in the tolerances /
##' heterogeneities be reduced via scaling by Hill's N2?
##' @param ... arguments passed to other methods
##' @return matrix of tolerances/heterogeneities with some additional
##'   attributes.
##' @author Gavin Simpson \email{gavin.simpson AT ucl.ac.uk}
##' @examples
##' data(dune)
##' data(dune.env)
##' mod <- cca(dune ~ ., data = dune.env)
##' tolerance.cca(mod)
##'
tolerance.cca <- function(x, choices = 1:2,
                          which = c("species","sites"),
                          scaling = 2, useN2 = FALSE, ...) {
    if(inherits(x, "rda"))
        stop("tolerances only available for unimodal ordinations")
    if(missing(which))
        which <- "species"
    ## reconstruct species/response matrix Y - up to machine precision!
    partialFit <- ifelse(is.null(x$pCCA$Fit), 0, x$pCCA$Fit)
    Y <- ((partialFit + x$CCA$Xbar) * sqrt(x$rowsum %o% x$colsum) +
          x$rowsum %o% x$colsum) * x$grand.total
    which <- match.arg(which)
    siteScrTypes <- if(is.null(x$CCA)){ "sites" } else {"lc"}
    scrs <- scores(x, display = c(siteScrTypes,"species"),
                   choices = choices, scaling = scaling)
    ## compute N2 if useN2 == TRUE & only if
    doN2 <- isTRUE(useN2) && ((which == "species" && abs(scaling) == 2) ||
                              (which == "sites" && abs(scaling) == 1))

    ## this gives the x_i - u_k on axis j
    ## outer(scrs$sites, scrs$species, "-")[,2,,j]
    siteScrs <- which(names(scrs) %in% c("sites","constraints"))
    xiuk <- outer(scrs[[siteScrs]], scrs$species, "-")
    if(isTRUE(all.equal(which, "sites"))) {
        ## need to permute the array as rowSums has different idea of what rows
        ## are that doesn't correspond to colSums. So flip dimensions 1 and 2
        ## with aperm and use colSums.
        res <- sqrt(sweep(colSums(aperm(sweep(xiuk[ , 2, , choices]^2, c(1:2),
                                              data.matrix(Y), "*"),
                                        c(2,1,3))),
                          1, rowSums(Y), "/"))
        if(doN2) {
            tot <- rowSums(Y)
            y <- sweep(Y, 1, tot, "/")^2
            N2 <- 1 / rowSums(y, na.rm = TRUE) ## 1/H
            res <- sweep(res, 1, sqrt(1 - (1/N2)), "/")
        }
    } else {
        res <- sqrt(sweep(colSums(sweep(xiuk[ , 2, , choices]^2, c(1:2),
                                        data.matrix(Y), "*")),
                          1, colSums(Y), "/"))
        if(doN2) {
            tot <- colSums(Y)
            y <- sweep(Y, 2, tot, "/")^2
            N2 <- 1 / colSums(y, na.rm = TRUE) ## 1/H
            res <- sweep(res, 1, sqrt(1 - (1/N2)), "/")
        }
    }
    class(res) <- c("tolerance.cca","tolerance","matrix")
    attr(res, "which") <- which
    attr(res, "scaling") <- scaling
    attr(res, "N2") <- NULL
    if(doN2)
        attr(res, "N2") <- N2
    attr(res, "model") <- deparse(substitute(mod))
    return(res)
}

`print.tolerance.cca` <- function(x, ...) {
    cat("\n")
    msg <- ifelse(attr(x, "which") == "species", "Species Tolerances",
                  "Sample Heterogeneities")
    writeLines(strwrap(msg, prefix = "\t"), sep = "\n\n")
    msg <- paste("Scaling:", attr(x, "scaling"))
    writeLines(strwrap(msg), sep = "\n\n")
    attr(x, "model") <- attr(x, "scaling") <- attr(x, "which") <- NULL
    print(unclass(x), ...)
    cat("\n")
}

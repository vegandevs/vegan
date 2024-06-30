`inertcomp` <-
    function (object, display = c("species", "sites"),
              unity = FALSE, proportional = FALSE)
{
    display <- match.arg(display)
    ## unity and proportional are conflicting arguments
    if (unity && proportional)
        stop("arguments 'unity' and 'proportional' cannot be both TRUE")
    if (!inherits(object, "cca"))
        stop("can be used only with objects inheriting from 'cca'")
    if (inherits(object, c("capscale", "dbrda")) && display == "species")
        stop(gettextf("cannot analyse species with '%s'", object$method))
    if (inherits(object, "dbrda")) {
        display <- "dbrda"
    }
    ## function to get the eigenvalues
    getComps <- function(x, display) {
        if(!is.null(x))
            switch(display,
                   "species" = colSums(x^2),
                   "sites" = rowSums(x^2),
                   "dbrda" = diag(x)
                   )
    }
    pCCA <- ordiYbar(object, "pCCA")
    CCA <- ordiYbar(object, "CCA")
    CA <- ordiYbar(object, "CA")
    tot <- ordiYbar(object, "initial")
    out <- cbind("pCCA" = getComps(pCCA, display),
                 "CCA" = getComps(CCA, display),
                 "CA" = getComps(CA, display))
    if (unity) ## each column sums to 1
        out <- sweep(out, 2, colSums(out), "/")
    if (proportional)
        out <- sweep(out, 1, rowSums(out), "/")
    out
}

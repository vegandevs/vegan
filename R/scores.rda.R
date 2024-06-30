### extract scores from rda, capscale and dbrda results. The two
### latter can have special features which are commented below. cca
### results are handled by scores.cca.
`scores.rda` <-
    function (x, choices = c(1, 2), display = "all", scaling = "species",
              const, correlation = FALSE, tidy = FALSE, droplist = TRUE,
              ...)
{
    ## Check the na.action, and pad the result with NA or WA if class
    ## "exclude"
    if (!is.null(x$na.action) && inherits(x$na.action, "exclude"))
        x <- ordiNApredict(x$na.action, x)
    tabula <- c("species", "sites", "constraints", "biplot",
                "regression", "centroids")
    names(tabula) <- c("sp", "wa", "lc", "bp", "reg", "cn")
    if (is.null(x$CCA))
        tabula <- tabula[1:2]
    display <- match.arg(display, c("sites", "species", "wa",
                                    "lc", "bp", "cn", "reg", "all"),
                         several.ok = TRUE)
    ## set "all" for tidy scores
    if (tidy)
        display <- "all"
    if("sites" %in% display)
      display[display == "sites"] <- "wa"
    if("species" %in% display)
      display[display == "species"] <- "sp"
    if ("all" %in% display)
        display <- names(tabula)
    take <- tabula[display]
    sumev <- x$tot.chi
    ## dbrda can have negative eigenvalues, but have scores only for
    ## positive
    eigval <- eigenvals(x)
    if (inherits(x, "dbrda") && any(eigval < 0))
        eigval <- eigval[eigval > 0]
    nr <- if (is.null(x$CCA))
        nrow(x$CA$u)
    else
        nrow(x$CCA$u)
    ## const multiplier of scores
    if (missing(const))
        const <- sqrt(sqrt((nr-1) * sumev))
    ## canoco 3 compatibility -- canoco 4 is incompatible
    ##else if (pmatch(const, "canoco")) {
    ##    const <- (sqrt(nr-1), sqrt(nr))
    ##}
    ##
    ## const[1] for species, const[2] for sites and friends
    if (length(const) == 1) {
        const <- c(const, const)
    }
    ## in dbrda we only have scores for positive eigenvalues
    if (inherits(x, "dbrda"))
        rnk <- x$CCA$poseig
    else
        rnk <- x$CCA$rank
    if (is.null(rnk))
        rnk <- 0
    choices <- choices[choices %in% seq_len(rnk + x$CA$rank)]
    slam <- sqrt(eigval[choices]/sumev)
    sol <- list()
    ## process scaling; numeric scaling will just be returned as is
    scaling <- scalingType(scaling = scaling, correlation = correlation)
    if ("species" %in% take) {
        v <- cbind(x$CCA$v, x$CA$v)[, choices, drop=FALSE]
        if (scaling) {
            scal <- list(1, slam, sqrt(slam))[[abs(scaling)]]
            v <- sweep(v, 2, scal, "*")
            if (scaling < 0) {
                v <- sweep(v, 1, x$colsum, "/")
                v <- v * sqrt(sumev / (nr - 1))
            }
            v <- const[1] * v
        }
        if (nrow(v) > 0)
            sol$species <- v
        else
            sol$species <- NULL
    }
    if ("sites" %in% take) {
        wa <- cbind(x$CCA$wa, x$CA$u)[, choices, drop=FALSE]
        if (scaling) {
            scal <- list(slam, 1, sqrt(slam))[[abs(scaling)]]
            wa <- sweep(wa, 2, scal, "*")
            wa <- const[2] * wa
        }
        sol$sites <- wa
    }
    if ("constraints" %in% take) {
        u <- cbind(x$CCA$u, x$CA$u)[, choices, drop=FALSE]
        if (scaling) {
            scal <- list(slam, 1, sqrt(slam))[[abs(scaling)]]
            u <- sweep(u, 2, scal, "*")
            u <- const[2] * u
        }
        sol$constraints <- u
    }
    if ("biplot" %in% take && !is.null(x$CCA$biplot)) {
        b <- matrix(0, nrow(x$CCA$biplot), length(choices))
        b[, choices <= rnk] <- x$CCA$biplot[, choices[choices <= rnk]]
        colnames(b) <- c(colnames(x$CCA$u), colnames(x$CA$u))[choices]
        rownames(b) <- rownames(x$CCA$biplot)
        if (scaling) {
            scal <- list(slam, 1, sqrt(slam))[[abs(scaling)]]
            b <- sweep(b, 2, scal, "*")
        }
        sol$biplot <- b
    }
    if ("regression" %in% take) {
        b <- coef(x, norm = TRUE)
        reg <- matrix(0, nrow(b), length(choices))
        reg[, choices <= rnk] <- b[, choices[choices <= rnk]]
        dimnames(reg) <- list(rownames(b),
                              c(colnames(x$CCA$u), colnames(x$CA$u))[choices])
        if (scaling) {
            scal <- list(slam, 1, sqrt(slam))[[abs(scaling)]]
            reg <- sweep(reg, 2, scal, "*")
        }
        sol$regression <- reg
    }
    if ("centroids" %in% take) {
        if (is.null(x$CCA$centroids))
            sol$centroids <- NULL
        else {
            cn <- matrix(0, nrow(x$CCA$centroids), length(choices))
            cn[, choices <= rnk] <- x$CCA$centroids[, choices[choices <= rnk]]
            colnames(cn) <- c(colnames(x$CCA$u), colnames(x$CA$u))[choices]
            rownames(cn) <- rownames(x$CCA$centroids)
            if (scaling) {
                scal <- list(slam, 1, sqrt(slam))[[abs(scaling)]]
                cn <- sweep(cn, 2, scal, "*")
                cn <- const[2] * cn
            }
            sol$centroids <- cn
        }
    }
    ## Take care that scores have names
    if (length(sol)) {
        for (i in seq_along(sol)) {
            if (is.matrix(sol[[i]]))
                rownames(sol[[i]]) <-
                    rownames(sol[[i]], do.NULL = FALSE,
                             prefix = substr(names(sol)[i], 1, 3))
        }
    }
    ## tidy scores
    if (tidy) {
        if (length(sol) == 0) # no scores
            return(NULL)
        ## re-group biplot arrays duplicating factor centroids
        if (!is.null(sol$biplot) && !is.null(sol$centroids)) {
            dup <- rownames(sol$biplot) %in% rownames(sol$centroids)
            if (any(dup)) {
                sol$factorbiplot <- sol$biplot[dup,, drop=FALSE]
                sol$biplot <- sol$biplot[!dup,, drop=FALSE]
            }
        }
        group <- sapply(sol, nrow)
        group <- rep(names(group), group)
        sol <- do.call(rbind, sol)
        label <- rownames(sol)
        sol <- as.data.frame(sol)
        sol$score <- as.factor(group)
        sol$label <- label
    }
    ## collapse const if both items identical
    if (identical(const[1], const[2]))
        const <- const[1]
    ## return NULL for list(), matrix for single scores, and a list
    ## for several scores
    if (droplist)
        sol <- switch(min(2, length(sol)), sol[[1]], sol)
    if (!is.null(sol))
        attr(sol, "const") <- const
    sol
}

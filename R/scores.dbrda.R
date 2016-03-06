### dbrda() has no 'sp' or 'wa' scores. We siletly ignore 'sp' and
### regard 'wa' or 'sites' as synonym of 'lc'. Some eigenvalues may be
### (and usually) are negative, and we also skip those axes.
`scores.dbrda` <-
    function (x, choices = c(1, 2), display = c("lc", "cn"),
              scaling = "sites", const, ...)
{
    ## 'const' can be a vector of length 2 in rda, but we accept only
    ## one value. In rda, the second item is for sites & friends: take
    ## only that if two given, and stop with error if there are more
    ## than two items of const
    if (missing(const))
        const <- 1
    if (length(const) > 1) {
        if (length(const) == 2)
            const <- const[2]
        else
            stop("'const' vector too long")
    }
    ## Check the na.action, and pad the result with NA or WA if class
    ## "exclude"
    if (!is.null(x$na.action) && inherits(x$na.action, "exclude"))
        x <- ordiNApredict(x$na.action, x)
    tabula <- c("sites", "constraints", "biplot", 
                "centroids")
    names(tabula) <- c("wa", "lc", "bp", "cn")
    if (is.null(x$CCA)) 
        tabula <- tabula[1:2]
    display <- match.arg(display, c("sites", "species", "wa",
                                    "lc", "bp", "cn"),
                         several.ok = TRUE)
    if("sites" %in% display)
        display[display == "sites"] <- "lc"
    if ("wa" %in% display)
        display[display == "wa"] <- "lc"
    take <- tabula[display]
    sumev <- x$tot.chi
    slam <- c(x$CCA$eig, x$CA$eig)
    slam <- slam[slam > 0]
    slam <- sqrt(slam)[choices]
    nr <- nobs(x)
    ## canoco 3 compatibility -- canoco 4 is incompatible
    ##else if (pmatch(const, "canoco")) {
    ##    const <- (sqrt(nr-1), sqrt(nr))
    ##}
    ##
    ## const[1] for species, const[2] for sites and friends
    rnk <- x$CCA$rank
    sol <- list()
    ## check scaling for character & process it if so. Numeric scaling
    ## is similar as in rda, but negative scaling do not make sense,
    ## and we only take positive values. Obviously scaling = 3 neither
    ## makes sense, but we let it be as a user problem.
    if (is.character(scaling)) {
        scaling <- abs(scalingType(scaling = scaling))
    }
    if ("constraints" %in% take) {
        u <- cbind(x$CCA$u, x$CA$u)
        choices <- choices[choices <= ncol(u)]
        u <- u[, choices, drop=FALSE]
        if (scaling) {
            scal <- list(slam, 1, sqrt(slam))[[abs(scaling)]]
            u <- sweep(u, 2, scal, "*")
        }
        u <- u * const
        sol$constraints <- u
    }
    if ("biplot" %in% take && !is.null(x$CCA$biplot)) {
        b <- matrix(0, nrow(x$CCA$biplot), length(choices))
        b[, choices <= rnk] <- x$CCA$biplot[, choices[choices <= rnk]]
        colnames(b) <- c(colnames(x$CCA$u), colnames(x$CA$u))[choices]
        rownames(b) <- rownames(x$CCA$biplot)
        sol$biplot <- b
    }
    if ("centroids" %in% take) {
        if (is.null(x$CCA$centroids))
            sol$centroids <- NA
        else {
            cn <- matrix(0, nrow(x$CCA$centroids), length(choices))
            cn[, choices <= rnk] <- x$CCA$centroids[, choices[choices <= rnk]]
            colnames(cn) <- c(colnames(x$CCA$u), colnames(x$CA$u))[choices]
            rownames(cn) <- rownames(x$CCA$centroids)
            if (scaling) {
                scal <- list(slam, 1, sqrt(slam))[[abs(scaling)]]
                cn <- sweep(cn, 2, scal, "*")
            }
            cn <- cn * const
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
    ## Only one type of scores: return a matrix instead of a list
    if (length(sol) == 1)
        sol <- sol[[1]]
    attr(sol, "const") <- const
    sol
}

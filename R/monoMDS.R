`monoMDS` <-
    function(dist, y, k = 2,
             model = c("global", "local", "linear", "hybrid"),
             threshold = 0.8, maxit = 200, weakties = TRUE, stress = 1,
             scaling = TRUE, pc = TRUE, smin = 1e-4, sfgrmin = 1e-7,
             sratmax=0.99999, ...) 
{
    ## Check that 'dist' are distances or a symmetric square matrix
    if (!(inherits(dist, "dist") ||
              ((is.matrix(dist) || is.data.frame(dist)) &&
                   isSymmetric(unname(as.matrix(dist))))))
        stop("'dist' must be a distance object (class \"dist\") or a symmetric square matrix")
    if (any(dist < -sqrt(.Machine$double.eps), na.rm = TRUE))
        warning("some dissimilarities are negative -- is this intentional?")
    ## match.arg
    model <- match.arg(model)
    ## save 'dist' attributes to display in print()
    distmethod <- attr(dist, "method")
    if (is.null(distmethod))
        distmethod <- "unknown"
    distcall <-attr(dist, "call")
    if (!is.null(distcall))
        distcall <- deparse(distcall)
    ## dist to mat
    mat <- as.matrix(dist)
    nm <- rownames(mat)
    if (model %in% c("global", "linear")) {
        ## global NMDS: lower triangle
        dist <- mat[lower.tri(mat)]
        iidx <- row(mat)[lower.tri(mat)]
        jidx <- col(mat)[lower.tri(mat)]
        ## Remove missing values
        if (any(nas <- is.na(dist))) {
            dist <- dist[!nas]
            iidx <- iidx[!nas]
            jidx <- jidx[!nas]
        }
        ## non-metric/metric: Fortran parameter 'iregn'
        if (model == "global")
            iregn <- 1
        else
            iregn <- 2
        ngrp <- 1
        nobj <- nrow(mat)
        istart <- 1
    } else if (model == "local") {
        ## local NMDS: whole matrix without the diagonal, and rows in
        ## a row (hence transpose)
        mat <- t(mat)
        ## Get missing values
        nas <- is.na(mat)
        ## groups by rows, except missing values
        rs <- rowSums(!nas) - 1
        istart <- cumsum(rs)
        istart <- c(0, istart[-length(istart)]) + 1
        ## Full matrix expect the diagonal
        dist <- mat[col(mat) != row(mat)]
        iidx <- col(mat)[col(mat) != row(mat)]  # transpose!
        jidx <- row(mat)[col(mat) != row(mat)]
        ## Remove missing values
        if (any(nas)) {
            nas <- nas[col(mat) != row(mat)]
            dist <- dist[!nas]
            iidx <- iidx[!nas]
            jidx <- jidx[!nas]
        }
        iregn <- 1
        nobj <- nrow(mat)
        ngrp <- nobj
    } else if (model == "hybrid") {
        ## Hybrid NMDS: two lower triangles, first a complete one,
        ## then those with dissimilarities below the threshold
        dist <- mat[lower.tri(mat)]
        iidx <- row(mat)[lower.tri(mat)]
        jidx <- col(mat)[lower.tri(mat)]
        ## Missing values
        if (any(nas <- is.na(dist))) {
            dist <- dist[!nas]
            iidx <- iidx[!nas]
            jidx <- jidx[!nas]
        }
        ## second group: dissimilarities below threshold
        ngrp <- 2
        istart <- c(1, length(dist) + 1)
        take <- dist < threshold
        dist <- c(dist, dist[take])
        iidx <- c(iidx, iidx[take])
        jidx <- c(jidx, jidx[take])
        iregn <- 3
        nobj <- nrow(mat)
    }
    ## ndis: number of >0 dissimilarities (distinct points)
    ndis <- length(dist)
    ## starting configuration
    if (missing(y)) {
        y <- matrix(runif(nobj*k, -1, 1), nobj, k)
        ## centre
        y <- sweep(y, 2, colMeans(y), "-")
    }
    ## y to vector
    y <- as.vector(as.matrix(y))
    ## translate R args to Fortran call
    if (weakties)
        ities <- 1
    else
        ities <- 2
    ## Fortran call
    sol <- .Fortran("monoMDS", nobj = as.integer(nobj), nfix=as.integer(0),
                 ndim = as.integer(k), ndis = as.integer(ndis),
                 ngrp = as.integer(ngrp), diss = as.double(dist),
                 iidx = as.integer(iidx), jidx = as.integer(jidx),
                 xinit = as.double(y), istart = as.integer(istart),
                 isform = as.integer(stress), ities = as.integer(ities),
                 iregn = as.integer(iregn), iscal = as.integer(scaling),
                 maxits = as.integer(maxit),
                 sratmx = as.double(sratmax), strmin = as.double(smin),
                 sfgrmn = as.double(sfgrmin), dist = double(ndis),
                 dhat = double(ndis), points = double(k*nobj),
                 stress = double(1), grstress = double(ngrp),
                 iters = integer(1), icause = integer(1),
                 PACKAGE = "vegan")
    sol$call <- match.call()
    sol$model <- model
    sol$points <- matrix(sol$points, nobj, k)
    if (pc)
        sol$points <- prcomp(sol$points)$x
    attr(sol$points, "pc") <- pc
    rownames(sol$points) <- nm
    colnames(sol$points) <- paste("MDS", 1:k, sep="")
    ## save info on dissimilarities
    sol$distmethod <- distmethod
    sol$distcall <- distcall
    class(sol) <- "monoMDS"
    sol
}


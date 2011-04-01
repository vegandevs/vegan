monoMDS <-
    function(dist, y, k = 2,
             model = c("global", "local", "hybrid"), threshold = 0.8,
             maxit = 200, tol = 0.0001, ...) 
{
    model <- match.arg(model)
    if (model == "global") {
        ## global NMDS: lower triangle
        mat <- as.matrix(dist)
        dist <- mat[lower.tri(mat)]
        iidx <- row(mat)[lower.tri(mat)]
        jidx <- col(mat)[lower.tri(mat)]
        iregn <- 1
        ngrp <- 1
        nobj <- nrow(mat)
        istart <- 1
    } else if (model == "local") {
        ## local NMDS: whole matrix without the diagonal, and rows in
        ## a row (hence transpose)
        mat <- t(as.matrix(dist))
        dist <- mat[col(mat) != row(mat)]
        iidx <- col(mat)[col(mat) != row(mat)]  # transpose!
        jidx <- row(mat)[col(mat) != row(mat)]
        iregn <- 1
        nobj <- nrow(mat)
        ngrp <- nobj
        istart <- seq(1, length(dist), by = (nobj-1))
    } else if (model == "hybrid") {
        ## Hybrid NMDS: two lower triangles, first a complete one,
        ## then those with dissimilarities below the threshold
        mat <- as.matrix(dist)
        dist <- mat[lower.tri(mat)]
        iidx <- row(mat)[lower.tri(mat)]
        jidx <- col(mat)[lower.tri(mat)]
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
    ## ndis: number dissimilarities
    ndis <- length(dist)
    ## starting configuration
    if (missing(y)) {
        y <- matrix(runif(nobj*k, -1, 1), nobj, k)
        ## centre
        y <- sweep(y, 2, colMeans(y), "-")
    }
    ## y to vector
    y <- as.vector(as.matrix(y))
    ## Fortran call
    sol <- .Fortran("monoMDS", nobj = as.integer(nobj), nfix=as.integer(0),
                 ndim = as.integer(k), ndis = as.integer(ndis),
                 ngrp = as.integer(ngrp), diss = as.double(dist),
                 iidx = as.integer(iidx), jidx = as.integer(jidx),
                 xinit = as.double(y), istart = as.integer(istart),
                 isform = as.integer(1), ities = as.integer(1),
                 iregn = as.integer(iregn), iscal = as.integer(1),
                 maxits = as.integer(maxit),
                 sratmx = as.double(0.99999), strmin = as.double(tol),
                 sfgrmn = as.double(1e-7), dist = double(ndis),
                 dhat = double(ndis), points = double(k*nobj),
                 stress = double(1), iters = integer(1),
                 icause = integer(1), PACKAGE = "vegan")
    sol$points <- matrix(sol$points, nobj, k)
    class(sol) <- "monoMDS"
    sol
}


`decorana` <-
    function (veg, iweigh = 0, iresc = 4, ira = 0, mk = 26, short = 0,
              before = NULL, after = NULL)
{
    ## constants
    Const2 <- 5
    Const3 <- 1e-11
    ZEROEIG <- 1e-7 # same limit as in the C function do_decorana
    ## data
    veg <- as.matrix(veg)
    if (!is.numeric(veg))
        stop("data 'veg' must be numeric (not factors or characters)")
    if (any(veg < 0))
        stop("'decorana' cannot handle negative data entries")
    ## optional data transformation
    if (!is.null(before)) {
        veg <- beforeafter(veg, before, after)
    }
    if (iweigh) {
        veg <- downweight(veg, Const2)
    }
    v <- attr(veg, "v")
    v.fraction <- attr(veg, "fraction")
    ## marginal sums after optional data transformations
    aidot <- rowSums(veg)
    if (any(aidot <= 0))
        stop("all row sums must be >0 in the community matrix: remove empty sites")
    adotj <- colSums(veg)
    if (any(adotj <= 0))
        warning("some species were removed because they were missing in the data")
    adotj[adotj < Const3] <- Const3
    ## check arguments
    if (mk < 10)
        mk <- 10
    if (mk > 46)
        mk <- 46
    if (ira)
        iresc <- 0
    ## Start analysis
    CA <- .Call(do_decorana, veg, ira, iresc, short, mk, as.double(aidot),
                as.double(adotj))
    if (ira)
        dnames <- paste("RA", 1:4, sep = "")
    else dnames <- paste("DCA", 1:4, sep = "")
    dimnames(CA$rproj) <- list(rownames(veg), dnames)
    dimnames(CA$cproj) <- list(colnames(veg), dnames)
    names(CA$evals) <- dnames
    origin <- apply(CA$rproj, 2, weighted.mean, aidot)
    vegChi <- initCA(veg) # needed for eigenvalues & their sum
    totchi <- sum(vegChi^2)
    if (ira) {
        evals.decorana <- NULL
        evals.ortho <- NULL
    }
    else {
        evals.decorana <- CA$evals
        if (any(ze <- evals.decorana <= 0))
            CA$evals[ze] <- 0
        ## centred and weighted scores
        x0 <- scale(CA$rproj, center = origin, scale = FALSE)
        x0 <- sqrt(aidot/sum(aidot)) * x0
        y0 <- scale(CA$cproj, center = origin, scale = FALSE)
        y0 <- sqrt(adotj/sum(adotj)) * y0
        ## eigenvalue: shrinking of scores y0 --WA--> x0
        evals <- colSums(x0^2)/colSums(y0^2)
        evals[evals < ZEROEIG | !is.finite(evals)] <- 0
        CA$evals <- evals
        ## decorana finds row scores from species scores, and for
        ## additive eigenvalues we need orthogonalized species
        ## scores. Q of QR decomposition will be orthonormal and if we
        ## use it for calculating row scores, these directly give
        ## additive eigenvalues.
        qy <- qr.Q(qr(y0))
        evals.ortho <- numeric(4) # qy can have < 4 columns
        evals.ortho[seq_len(ncol(qy))] <- colSums(crossprod(t(vegChi), qy)^2)
        evals.ortho[evals.ortho < ZEROEIG | !is.finite(evals.ortho)] <- 0
        names(evals.ortho) <- names(evals.decorana)
    }
    additems <- list(totchi = totchi, evals.ortho = evals.ortho,
                     evals.decorana = evals.decorana, origin = origin,
                     v = v, fraction = v.fraction, iweigh = iweigh,
                     before = before, after = after,
                     call = match.call())
    CA <- c(CA, additems)
    class(CA) <- "decorana" # c() strips class
    CA
}

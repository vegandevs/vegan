`decorana2` <-
    function (veg, iweigh = 0, iresc = 4, ira = 0, mk = 26, short = 0,
              before = NULL, after = NULL)
{
    Const1 <- 1e-10
    Const2 <- 5
    Const3 <- 1e-11
    veg <- as.matrix(veg)
    if (any(rowSums(veg) <= 0))
        stop("all row sums must be >0 in the community matrix: remove empty sites")
    if (any(veg < 0))
        stop("'decorana' cannot handle negative data entries")
    if (any(colSums(veg) <= 0))
        warning("some species were removed because they were missing in the data")
    if (mk < 10)
        mk <- 10
    if (mk > 46)
        mk <- 46
    if (ira)
        iresc <- 0
    if (!is.null(before)) {
        if (is.unsorted(before))
            stop("'before' must be sorted")
        if (length(before) != length(after))
            stop("'before' and 'after' must have same lengths")
        for (i in 1:nr) {
            tmp <- veg[i, ] > 0
            veg[i, tmp] <- approx(before, after, veg[i, tmp],
                                  rule = 2)$y
        }
    }
    if (iweigh) {
        veg <- downweight(veg, Const2)
    }
    v <- attr(veg, "v")
    v.fraction <- attr(veg, "fraction")
    adotj <- colSums(veg)
    adotj[adotj < Const3] <- Const3
    aidot <- rowSums(veg)
    CA <- .Call("do_decorana", veg, ira, iresc, short, mk, as.double(aidot),
                as.double(adotj), PACKAGE = "vegan")
    if (ira)
        dnames <- paste("RA", 1:4, sep = "")
    else dnames <- paste("DCA", 1:4, sep = "")
    rownames(CA$rproj) <- rownames(veg)
    colnames(CA$rproj) <- dnames
    rownames(CA$cproj) <- colnames(veg)
    colnames(CA$cproj) <- dnames
    names(CA$evals) <- dnames
    origin <- apply(CA$rproj, 2, weighted.mean, aidot)
    if (ira) {
        evals.decorana <- NULL
    }
    else {
        evals.decorana <- CA$evals
        var.r <- diag(cov.wt(CA$rproj, aidot, method = "ML")$cov)
        var.c <- diag(cov.wt(CA$cproj, adotj, method = "ML")$cov)
        CA$evals <- var.r/var.c
        if (any(ze <- evals.decorana <= 0))
            CA$evals[ze] <- 0
    }
    CA <- list(rproj = CA$rproj, cproj = CA$cproj, evals = CA$evals,
               evals.decorana = evals.decorana,
               origin = origin, v = v, fraction = v.fraction, adotj = CA$adotj,
               aidot = CA$aidot, iweigh = iweigh, iresc = CA$iresc, ira = CA$ira,
               mk = CA$mk, short = CA$short, before = before, after = after,
               call = match.call())
    class(CA) <- "decorana"
    CA
}

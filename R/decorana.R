`decorana` <-
    function (veg, iweigh = 0, iresc = 4, ira = 0, mk = 26, short = 0,
              before = NULL, after = NULL)
{
    Const1 <- 1e-10
    Const2 <- 5
    Const3 <- 1e-11
    veg <- as.matrix(veg)
    aidot <- rowSums(veg)
    if (any(aidot <= 0))
        stop("all row sums must be >0 in the community matrix: remove empty sites")
    if (any(veg < 0))
        stop("'decorana' cannot handle negative data entries")
    adotj <- colSums(veg)
    if (any(adotj <= 0))
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
        for (i in seq_len(nrow(veg))) {
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
    adotj[adotj < Const3] <- Const3
    CA <- .Call(do_decorana, veg, ira, iresc, short, mk, as.double(aidot),
                as.double(adotj))
    if (ira)
        dnames <- paste("RA", 1:4, sep = "")
    else dnames <- paste("DCA", 1:4, sep = "")
    dimnames(CA$rproj) <- list(rownames(veg), dnames)
    dimnames(CA$cproj) <- list(colnames(veg), dnames)
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
    additems <- list(evals.decorana = evals.decorana, origin = origin,
                     v = v, fraction = v.fraction, iweigh = iweigh,
                     before = before, after = after,
                     call = match.call())
    CA <- c(CA, additems)
    class(CA) <- "decorana" # c() strips class
    CA
}

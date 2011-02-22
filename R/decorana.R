"decorana" <-
    function (veg, iweigh = 0, iresc = 4, ira = 0, mk = 26, short = 0, 
              before = NULL, after = NULL) 
{
    Const1 <- 1e-10
    Const2 <- 5
    Const3 <- 1e-11
    ZEROEIG <- 1e-7 # consider as zero eigenvalue
    veg <- as.matrix(veg)
    if (any(rowSums(veg) <= 0)) 
        stop("All row sums must be >0 in the community matrix: remove empty sites.")
    if (any(veg < 0))
        stop("'decorana' cannot handle negative data entries")
    if (any(colSums(veg) <= 0)) 
        warning("Some species were removed because they were missing in the data.")
    nr <- nrow(veg)
    nc <- ncol(veg)
    mk <- mk + 4
    if (mk < 14) 
        mk <- 14
    if (mk > 50) 
        mk <- 50
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
    tot <- sum(adotj)
    yeig1 <- rep(1, nc)
    xeig1 <- rep(1, nr)
    eig <- 1
    nid <- sum(veg > 0)
    cep <- .C("data2hill", as.double(veg), mi = as.integer(nr), 
              n = as.integer(nc), nid = as.integer(nid), ibegin = integer(nr), 
              iend = integer(nr), idat = integer(nid), qidat = double(nid), 
              PACKAGE = "vegan")[c("mi", "n", "nid", "ibegin", "iend", 
              "idat", "qidat")]
    ix1 <- ix2 <- ix3 <- rep(0, cep$mi)
    s1 <- .Fortran("eigy", x = as.double(xeig1), y = as.double(yeig1), 
                   eig = double(1), neig = as.integer(0), ira = as.integer(ira), 
                   iresc = as.integer(iresc), short = as.double(short), 
                   mi = as.integer(cep$mi), mk = as.integer(mk), n = as.integer(cep$n), 
                   nid = as.integer(cep$ni), ibegin = as.integer(cep$ibegin), 
                   iend = as.integer(cep$iend), idat = as.integer(cep$idat), 
                   qidat = as.double(cep$qidat), y2 = double(cep$n), y3 = double(cep$n), 
                   y4 = double(cep$n), y5 = double(cep$n), xeig1 = as.double(xeig1), 
                   xeig2 = double(cep$mi), xeig3 = double(cep$mi), ix1 = as.integer(ix1), 
                   ix2 = as.integer(ix2), ix3 = as.integer(ix3), aidot = as.double(aidot), 
                   adotj = as.double(adotj), PACKAGE = "vegan")[c("x", "y", 
                                             "eig")]
    ## Eigenvalues can be computed zeros: zero the results exactly to
    ## avoid problems with rescaling and estimating final eigenvalues
    ## in arbitrary corner cases.
    if (s1$eig < ZEROEIG)
        s1$x[] <- s1$y[] <- s1$eig <- 0
    if (!ira) 
        ix1 <- .Fortran("cutup", x = as.double(s1$x), ix = as.integer(ix1), 
                        mi = as.integer(cep$mi), mk = as.integer(mk), PACKAGE = "vegan")$ix
    s2 <- .Fortran("eigy", x = as.double(xeig1), y = as.double(yeig1), 
                   eig = double(1), neig = as.integer(1), ira = as.integer(ira), 
                   iresc = as.integer(iresc), short = as.double(short), 
                   mi = as.integer(cep$mi), mk = as.integer(mk), n = as.integer(cep$n), 
                   nid = as.integer(cep$ni), ibegin = as.integer(cep$ibegin), 
                   iend = as.integer(cep$iend), idat = as.integer(cep$idat), 
                   qidat = as.double(cep$qidat), y2 = double(cep$n), y3 = double(cep$n), 
                   y4 = double(cep$n), y5 = double(cep$n), xeig1 = as.double(s1$x), 
                   xeig2 = double(cep$mi), xeig3 = double(cep$mi), ix1 = as.integer(ix1), 
                   ix2 = as.integer(ix2), ix3 = as.integer(ix3), aidot = as.double(aidot), 
                   adotj = as.double(adotj), PACKAGE = "vegan")[c("x", "y", 
                                             "eig")]
    if (s2$eig < ZEROEIG)
        s2$x[] <- s2$y[] <- s2$eig <- 0
    if (!ira) 
        ix2 <- .Fortran("cutup", x = as.double(s2$x), ix = as.integer(ix2), 
                        mi = as.integer(cep$mi), mk = as.integer(mk), PACKAGE = "vegan")$ix
    s3 <- .Fortran("eigy", x = as.double(xeig1), y = as.double(yeig1), 
                   eig = double(1), neig = as.integer(2), ira = as.integer(ira), 
                   iresc = as.integer(iresc), short = as.double(short), 
                   mi = as.integer(cep$mi), mk = as.integer(mk), n = as.integer(cep$n), 
                   nid = as.integer(cep$ni), ibegin = as.integer(cep$ibegin), 
                   iend = as.integer(cep$iend), idat = as.integer(cep$idat), 
                   qidat = as.double(cep$qidat), y2 = double(cep$n), y3 = double(cep$n), 
                   y4 = double(cep$n), y5 = double(cep$n), xeig1 = as.double(s1$x), 
                   xeig2 = as.double(s2$x), xeig3 = double(cep$mi), ix1 = as.integer(ix1), 
                   ix2 = as.integer(ix2), ix3 = as.integer(ix3), aidot = as.double(aidot), 
                   adotj = as.double(adotj), PACKAGE = "vegan")[c("x", "y", 
                                             "eig")]
    if (s3$eig < ZEROEIG)
        s3$x[] <- s3$y[] <- s3$eig <- 0
    if (!ira) 
        ix3 <- .Fortran("cutup", x = as.double(s3$x), ix = as.integer(ix3), 
                        mi = as.integer(cep$mi), mk = as.integer(mk), PACKAGE = "vegan")$ix
    s4 <- .Fortran("eigy", x = as.double(xeig1), y = as.double(yeig1), 
                   eig = double(1), neig = as.integer(3), ira = as.integer(ira), 
                   iresc = as.integer(iresc), short = as.double(short), 
                   mi = as.integer(cep$mi), mk = as.integer(mk), n = as.integer(cep$n), 
                   nid = as.integer(cep$ni), ibegin = as.integer(cep$ibegin), 
                   iend = as.integer(cep$iend), idat = as.integer(cep$idat), 
                   qidat = as.double(cep$qidat), y2 = double(cep$n), y3 = double(cep$n), 
                   y4 = double(cep$n), y5 = double(cep$n), xeig1 = as.double(s1$x), 
                   xeig2 = as.double(s2$x), xeig3 = as.double(s3$x), ix1 = as.integer(ix1), 
                   ix2 = as.integer(ix2), ix3 = as.integer(ix3), aidot = as.double(aidot), 
                   adotj = as.double(adotj), PACKAGE = "vegan")[c("x", "y", 
                                             "eig")]
    if (s4$eig < ZEROEIG)
        s4$x[] <- s4$y[] <- s4$eig <- 0
    if (s1$eig > ZEROEIG)
        s1$x <- .Fortran("yxmult", y = as.double(s1$y), x = as.double(s1$x), 
                     as.integer(cep$mi), as.integer(cep$n), as.integer(cep$nid), 
                     as.integer(cep$ibegin), as.integer(cep$iend), as.integer(cep$idat), 
                     as.double(cep$qidat), PACKAGE = "vegan")$x/aidot
    if (s2$eig > ZEROEIG)
        s2$x <- .Fortran("yxmult", y = as.double(s2$y), x = as.double(s2$x), 
                     as.integer(cep$mi), as.integer(cep$n), as.integer(cep$nid), 
                     as.integer(cep$ibegin), as.integer(cep$iend), as.integer(cep$idat), 
                     as.double(cep$qidat), PACKAGE = "vegan")$x/aidot
    if (s3$eig > ZEROEIG)
        s3$x <- .Fortran("yxmult", y = as.double(s3$y), x = as.double(s3$x), 
                     as.integer(cep$mi), as.integer(cep$n), as.integer(cep$nid), 
                     as.integer(cep$ibegin), as.integer(cep$iend), as.integer(cep$idat), 
                     as.double(cep$qidat), PACKAGE = "vegan")$x/aidot
    if (s4$eig > ZEROEIG)
        s4$x <- .Fortran("yxmult", y = as.double(s4$y), x = as.double(s4$x), 
                     as.integer(cep$mi), as.integer(cep$n), as.integer(cep$nid), 
                     as.integer(cep$ibegin), as.integer(cep$iend), as.integer(cep$idat), 
                     as.double(cep$qidat), PACKAGE = "vegan")$x/aidot
    rproj <- cbind(s1$x, s2$x, s3$x, s4$x)
    cproj <- cbind(s1$y, s2$y, s3$y, s4$y)
    evals <- c(s1$eig, s2$eig, s3$eig, s4$eig)
    if (ira) 
        dnames <- paste("RA", 1:4, sep = "")
    else dnames <- paste("DCA", 1:4, sep = "")
    rownames(rproj) <- rownames(veg)
    colnames(rproj) <- dnames
    rownames(cproj) <- colnames(veg)
    colnames(cproj) <- dnames
    names(evals) <- dnames
    origin <- apply(rproj, 2, weighted.mean, aidot)
    if (ira) {
        evals.decorana <- NULL
    }
    else {
        evals.decorana <- evals
        var.r <- cov.wt(rproj, aidot)
        var.r <- diag(var.r$cov) * (1 - sum(var.r$wt^2))
        var.c <- cov.wt(cproj, adotj)
        var.c <- diag(var.c$cov) * (1 - sum(var.c$wt^2))
        evals <- var.r/var.c
        if (any(ze <- evals.decorana < ZEROEIG))
            evals[ze] <- 0
    }
    CA <- list(rproj = rproj, cproj = cproj, evals = evals, evals.decorana = evals.decorana, 
               origin = origin, v = v, fraction = v.fraction, adotj = adotj, 
               aidot = aidot, iweigh = iweigh, iresc = iresc, ira = ira, 
               mk = mk - 4, short = short, before = before, after = after, 
               call = match.call())
    class(CA) <- "decorana"
    CA
}

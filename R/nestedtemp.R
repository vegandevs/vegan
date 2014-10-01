`nestedtemp` <-
    function(comm, ...)
{
    ## J Biogeogr 33, 924-935 (2006) says that Atmar & Patterson try
    ## to pack presences and absence to minimal matrix temperature,
    ## and the following routines try to reproduce the (partly verbal)
    ## description. Index s should pack ones, and index t should pack
    ## zeros, and the final ordering should be "a compromise".
    colpack <- function(x, rr)
    {
        ind <- matrix(rep(rr, ncol(x)), nrow=nrow(x))
        s <- -colSums((x*ind)^2)
        t <- -colSums((nrow(x) - (1-x)*ind + 1)^2)
        st <- rank(s+t, ties.method = "random")
        st
    }
    rowpack <- function(x, cr)
    {
        ind <- matrix(rep(cr, each=nrow(x)), nrow=nrow(x))
        s <- -rowSums((x*ind)^2)
        t <- -rowSums((ncol(x) - (1-x)*ind + 1)^2)
        st <- rank(s+t, ties.method = "random")
        st
    }
    comm <- ifelse(comm > 0, 1, 0)
    ## Start with columns, expect if nrow > ncol
    if (ncol(comm) >= nrow(comm)) {
        i <- rank(-rowSums(comm), ties.method = "average")
    } else {
        j <- rank(-colSums(comm), ties.method = "average")
        i <- rowpack(comm, j)
    }
    ## Improve eight times
    for (k in 1:8) {
        j <- colpack(comm, i)
        i <- rowpack(comm, j)
    }
    if (ncol(comm) < nrow(comm))
        j <- colpack(comm, i)
    comm <- comm[order(i), order(j)]
    r <- ppoints(nrow(comm), a=0.5)
    c <- ppoints(ncol(comm), a=0.5)
    dis <- matrix(rep(r, ncol(comm)), nrow=nrow(comm))
    totdis <- 1 - abs(outer(r, c, "-"))
    fill <- sum(comm)/prod(dim(comm))
    ## Fill line as defined in J Biogeogr by solving an integral of
    ## the fill function
    fillfun <- function(x, p) 1 - (1-(1-x)^p)^(1/p)
    intfun <- function(p, fill)
        integrate(fillfun, lower=0, upper=1, p=p)$value - fill
    ## 'p' will depend on 'fill', and fill = 0.0038 correspond to p =
    ## 20. Sometimes the fill is lower, and therefore we try() to see
    ## if we need to move the bracket up. We should need to do this
    ## very rarely.
    lo <- 0
    hi <- 20
    repeat{
        sol <- try(uniroot(intfun, c(lo,hi), fill=fill), silent = TRUE)
        if (inherits(sol, "try-error")) {
            if (hi > 640) # bail out
                stop(gettextf("matrix is too sparse, fill is %g"), fill)
            lo <- hi
            hi <- hi + hi
        } else
            break
    }
    p <- sol$root
    ## row coordinates of the fill line for all matrix entries
    out <- matrix(0, nrow=length(r), ncol=length(c))
    for (i in 1:length(r))
        for (j in 1:length(c)) {
            a <- c[j] - r[i]
            out[i,j] <- uniroot(function(x, ...) fillfun(x, p) - a -x,
                                c(0,1), p = p)$root
            }
    ## Filline
    x <- seq(0,1,len=51)
    xline <- fillfun(x, p)
    smo <- list(x = x, y = xline)
    u <- (dis - out)/totdis
    u[u < 0 & comm == 1] <- 0
    u[u > 0 & comm == 0] <- 0
    u <- u^2
    colnames(u) <- colnames(comm)
    rownames(u) <- rownames(comm)
    names(r) <- rownames(comm)
    names(c) <- colnames(comm)
    temp <- 100*sum(u)/prod(dim(comm))/0.04145
    out <- list(comm = comm, u = u, r = r, c = c, p = p,
                fill=fill,  statistic = temp, smooth=smo)
    names(out$statistic) <- "temperature"
    class(out) <- "nestedtemp"
    out
}


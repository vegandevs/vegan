## permatswap function
`permatswap` <-
function(m, method="quasiswap", reg=NULL, hab=NULL, mtype="count", times=100, burnin = 10000, thin = 1000)
{

## temporary internal function
## quasiswapcount is based on the idea of Carsten Dorman from function swap.web::bipartite
## this swapcount implementation roughly follow the original swapcount that was put in C
## new lines are commented
swapcount <-
function(m, thin = 1, mfill = 0)                                            # new arg
{
## internal, is the 2x2 matrix diagonal or anti-diagonal
## 'isDiag' is a utility function for 'swapcount' to find the largest
## value that can be swapped and whether in diagonal or antidiagonal
## way. The input is a 2x2 submatrix.
isDiag <- function(x) {
        x<- as.vector(x)
        X <- as.numeric(x>0)
        ## sX: number of non-zero cells
        sX <- sum(X)
        ## Smallest diagonal and antidiagonal element
        choose <- c(min(x[c(2,3)]), min(x[c(1,4)]))
        ## Either choose could be returned, but RNG is not needed,
        ## because submatrix already is in random order, and we always return choose[0]
        if (sX == 4) return(choose[1])
        if (identical(X, c(0,1,1,0)) || identical(X, c(0,1,1,1)) || identical(X, c(1,1,1,0)))
                return(choose[1])
        if (identical(X, c(1,0,0,1)) || identical(X, c(1,0,1,1)) || identical(X, c(1,1,0,1)))
                return(-choose[2])
        if (sX < 2 | identical(X, c(0,0,1,1)) || identical(X, c(1,1,0,0)) || 
            identical(X, c(0,1,0,1)) || identical(X, c(1,0,1,0)))
                return(0)
        } ## end of isDiag

    x <- as.matrix(m)
    n.col <- ncol(x)
    n.row <- nrow(x)
    changed <- 0
    if (mfill != 0) thin <- 1                                                                 # just to be sure
    while(changed < thin) {
        ran.row <- sample(n.row, 2)
        ran.col <- sample(n.col, 2)
        ## The largest value that can be swapped
        ev <- isDiag(x[ran.row, ran.col])
        if (ev != 0) {
            ## Check that the fill doesn't change
                                                                            # new condition: mfill == 0
            if (mfill == 0 && identical(sum(x[ran.row, ran.col] > 0), sum(x[ran.row, ran.col] + matrix(c(ev,-ev,-ev,ev), 2, 2) > 0)))
                {
                ## Swap
                x[ran.row, ran.col] <- x[ran.row, ran.col] + matrix(c(ev,-ev,-ev,ev), 2, 2)
                changed <- changed + 1
                } else {                                                                       # these few lies are also new
                                                                                               # logical: has it tecreased?
                decreased <- sum(x[ran.row, ran.col] > 0) > sum(x[ran.row, ran.col] + matrix(c(ev,-ev,-ev,ev), 2, 2) > 0)
                    # sum(x > 0) > mfill && decreased         DO
                    # sum(x > 0) > mfill && !decreased        DONT DO
                    # sum(x > 0) < mfill && !decreased         DO
                    # sum(x > 0) < mfill && decreased        DONT DO
                if ((sum(x > 0) > mfill && decreased) || (sum(x > 0) < mfill && !decreased))
                    # this result only depend on actual fill and ev effects
                    x[ran.row, ran.col] <- x[ran.row, ran.col] + matrix(c(ev,-ev,-ev,ev), 2, 2)
                if (sum(x > 0) == mfill)
                    changed <- changed + 1                         # changed changes only if target is met
                }                                                                   # new lines untill here
            }
        }
    return(x)
}



    if (!identical(all.equal(m, round(m)), TRUE))
       stop("function accepts only integers (counts)")
    mtype <- match.arg(mtype, c("prab", "count"))
    count <- mtype == "count"
    if (count) {
        method <- match.arg(method, c("swap", "quasiswap"))
    } else {method <- match.arg(method, c("swap", "quasiswap", "tswap", "backtracking"))}

    m <- as.matrix(m)
    n.row <- nrow(m)
    n.col <- ncol(m)
    if (mtype == "prab") m <- ifelse(m > 0, 1, 0)
    if (is.null(reg) && is.null(hab)) str <- as.factor(rep(1, n.row))
    if (!is.null(reg) && is.null(hab)) str <- as.factor(reg)
    if (is.null(reg) && !is.null(hab)) str <- as.factor(hab)
    if (!is.null(reg) && !is.null(hab)) str <- interaction(reg, hab, drop=TRUE)
    levels(str) <- 1:length(unique(str))
    str <- as.numeric(str)
    nstr <- length(unique(str))
    if (any(tapply(str,list(str),length) == 1))
        stop("strata should contain at least 2 observations")

    if (method != "quasiswap") {
        perm <- list()
        for (i in 1:times)
            perm[[i]] <- matrix(0, n.row, n.col)
    } else {
        perm <- r2dtable(times, apply(m, 1, sum), apply(m, 2, sum))
    }

    for (j in 1:nstr) {
        id <- which(str == j)
        temp <- m[id,]
        if (method != "quasiswap") {
            if (count)
                for (k in 1:burnin)
                    temp <- .C("swapcount", m = as.double(temp),
                            as.integer(n.row), as.integer(n.col),
                            as.integer(1), PACKAGE = "vegan")$m
            else
                for (k in 1:burnin)
                    temp <- commsimulator(temp, method=method)
            for (i in 1:times) {
                if (count)
                    perm[[i]][id,] <- .C("swapcount",
                                    m = as.double(temp),
                                    as.integer(n.row),
                                    as.integer(n.col),
                                    as.integer(thin),
                                    PACKAGE = "vegan")$m
	        else perm[[i]][id,] <- commsimulator(temp, method=method, thin=thin)
            temp <- perm[[i]][id,]
            } # for i end
        } else {
            for (i in 1:times) {
                if (count)
                    if (sum(perm[[i]][id,] > 0) != sum(m[id,] > 0))             ## if fills are equal, no need to do it (r2dtable)
                        perm[[i]][id,] <- swapcount(perm[[i]][id,], thin=1, mfill=sum(m[id,] > 0))   ## this should be replaced by .C()
                else perm[[i]][id,] <- commsimulator(temp, method=method)
            }
            thin <- burnin <- 0
        }
    } # for j end
    specs <- list(reg=reg, hab=hab, burnin=burnin, thin=thin)
    out <- list(call=match.call(), orig=m, perm=perm, specs=specs)
    attr(out, "mtype") <- mtype
    attr(out, "ptype") <- "swap"
    attr(out, "method") <- method
    attr(out, "fixedmar") <- "both"
    attr(out, "times") <- times
    attr(out, "shuffle") <- NA
    class(out) <- c("permat", "list")
    return(out)
}



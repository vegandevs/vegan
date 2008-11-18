"commsimulator" <-
function (x, method, thin = 1) 
{
    method <- match.arg(method, 
                        c("r0","r1","r2","r00","c0","swap", "tswap",
                          "backtrack", "quasiswap"))
    if (any(x > 1))
        x <- ifelse(x > 0, 1, 0)
    nr <- nrow(x)
    nc <- ncol(x)
    if (method %in% c("r0", "r1", "r2")) {
        rs <- rowSums(x)
        if (method == "r0")
            p <- rep(1, nc)
        else
            p <- colSums(x)
        if (method == "r2")
            p <- p*p
        out <- matrix(0, nrow=nr, ncol=nc)
        for (i in 1:nr)
            out[i,sample(nc, rs[i], prob=p)] <- 1 
    }
    else if (method == "r00") {
        out <- numeric(nr*nc)
        out[sample(length(out), sum(x))] <- 1
        dim(out) <- dim(x)
    }
    else if (method == "c0") {
        cs <- colSums(x)
        out <- matrix(0, nrow=nr, ncol=nc)
        for (j in 1:nc)
            out[sample(nr, cs[j]), j] <- 1
    } else if (method == "swap") {
        x <- as.matrix(x)
        out <- .C("swap", m = as.integer(x), as.integer(nrow(x)),
                  as.integer(ncol(x)), as.integer(thin),
                  PACKAGE = "vegan")$m
        dim(out) <- dim(x)
    } else if (method == "tswap") {
        x <- as.matrix(x)
        out <- .C("trialswap", m = as.integer(x), as.integer(nrow(x)),
                  as.integer(ncol(x)), as.integer(thin),
                  PACKAGE = "vegan")$m
        dim(out) <- dim(x)
    } else if (method == "quasiswap") {
        out <- r2dtable(1, rowSums(x), colSums(x))[[1]]
        out <- .C("quasiswap", m = as.integer(out), as.integer(nrow(x)),
                  as.integer(ncol(x)), PACKAGE = "vegan")$m
        dim(out) <- dim(x)
    }
    else if (method == "backtrack") {
        fill <- sum(x)
        rs <- rowSums(x) 
        cs <- colSums(x) 
        all <- matrix(1:(nr*nc), nrow=nr, ncol=nc)
        out <- matrix(0, nrow=nr, ncol=nc)
        free <- matrix(1:(nr*nc), nrow=nr)
        icount <- numeric(length(rs))
        jcount <- numeric(length(cs))
        ## Fill: ordering by cell probabilities
        prob <- outer(rs, cs, "*") 
        ij <- sample(free, prob=prob)
        i <- (ij - 1) %% nr + 1
        j <- (ij - 1) %/% nr + 1
        for (k in 1:length(ij)) {
            if (icount[i[k]] < rs[i[k]] && jcount[j[k]] < cs[j[k]]) {
            	out[ij[k]] <- 1
            	icount[i[k]] <- icount[i[k]] + 1
            	jcount[j[k]] <- jcount[j[k]] + 1
            }
        }
        ## "Backtrack": remove a random presence and fill with !present
        ndrop <- 1
        for (i in 1:10000) {
            oldout <- out
            oldn <- sum(out)
            drop <- sample(all[out==1], ndrop)
            out[drop] <- 0
            candi <- outer(rowSums(out) < rs, colSums(out) < cs, "&") & out == 0
            while (sum(candi) > 0) {
                if (sum(candi) > 1)
                    ij <- sample(all[candi], 1)
                else
                    ij <- all[candi]
                out[ij] <- 1
                candi <- outer(rowSums(out) < rs, colSums(out) < cs, "&") & out == 0
            }
            if (sum(out) >= fill) break
            if (oldn >= sum(out))
                ndrop <- min(ndrop + 1, 4)
            else
                ndrop <- 1
            if (oldn > sum(out))
                out <- oldout
        }
    }
    out
}


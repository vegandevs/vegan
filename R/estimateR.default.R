"estimateR.default" <-
    function (x, ...) 
{
    gradF <- function(a, i) {
        .expr4 <- sum(i * a)
        .expr7 <- 1 - a[1]/(1 - sum(i * a))
        .expr8 <- 1/.expr7
        .expr10 <- sum(a)
        .expr12 <- sum(i * (i - 1) * a)
        .expr13 <- sum(a) * sum(i * (i - 1) * a)
        .expr14 <- .expr7 * .expr4
        .expr15 <- .expr4 - 1
        .expr16 <- .expr14 * .expr15
        .expr18 <- .expr13/.expr16 - 1
        .expr20 <- sum(a) + a[1] * .expr18
        .expr23 <- (1 - sum(i * a))^2
        .expr25 <- 1/(1 - sum(i * a)) + a[1]/(1 - sum(i * a))^2
        .expr26 <- .expr7^2
        .expr35 <- .expr16^2
        Grad <- a[1] * i/(.expr23 * .expr26) * .expr20 + .expr8 * 
            (1 + a[1] * ((.expr12 + (.expr10 * i * (i - 1)))/.expr16 - 
                         .expr13 * ((.expr7 * i - (a[1] * i/.expr23) * 
                                     .expr4) * .expr15 + .expr14 * i)/.expr35))
        Grad[1] <- .expr25/.expr26 * .expr20 + .expr8 * (1 + 
                                                         (.expr18 + a[1] * (.expr12/.expr16 - .expr13 * ((.expr7 - 
                                                                                                          .expr25 * .expr4) * .expr15 + .expr14)/.expr35)))
        return(Grad)
    }
    if (!identical(all.equal(x, round(x)), TRUE)) 
        stop("function accepts only integers (counts)")
    freq <- x[x > 0]
    X <- x[x > 0]
    T.X <- table(X)
    S.obs <- length(X)
    S.rare <- sum(T.X[as.numeric(names(T.X)) <= 10])
    S.abund <- sum(T.X[as.numeric(names(T.X)) > 10])
    N.rare <- sum(X[X < 11])
    i <- 1:10
    COUNT <- function(i, counts) {
        length(counts[counts == i])
    }
    a <- sapply(i, COUNT, X)
    G <- a[1]/a[2]
    S.Chao1 <- S.obs + a[1] * (a[1] - 1) / (a[2] + 1)/ 2
    Deriv.Ch1 <- gradF(a, i)
    sd.Chao1 <- sqrt(a[2] * ((G^4)/4 + G^3 + (G^2)/2))
    C.ace <- 1 - a[1]/N.rare
    i <- 1:length(a)
    thing <- i * (i - 1) * a
    Gam <- sum(thing) * S.rare/(C.ace * N.rare * (N.rare - 1)) - 
        1
    S.ACE <- S.abund + S.rare/C.ace + max(Gam, 0) * a[1]/C.ace
    sd.ACE <- sqrt(sum(Deriv.Ch1 %*% t(Deriv.Ch1) * (diag(a) - 
                                                     a %*% t(a)/S.ACE)))
    out <- list(S.obs = S.obs, S.chao1 = S.Chao1, se.chao1 = sd.Chao1, 
                S.ACE = S.ACE, se.ACE = sd.ACE)
    out <- unlist(out)
    out
}

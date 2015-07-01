`estimateR.default` <-
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
        Grad
    }
    if (!identical(all.equal(x, round(x)), TRUE)) 
        stop("function accepts only integers (counts)")
    X <- x[x > 0]
    N <- sum(X)
    SSC <- 1 # (N-1)/N # do NOT use small-sample correction
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
    ## EstimateS uses basic Chao only if a[2] > 0, and switches to
    ## bias-corrected version only if a[2] == 0. However, we always
    ## use bias-corrected form. The switchin code is commented out so
    ## that it is easy to put back.

    ##if (a[2] > 0)
    ##    S.Chao1 <- S.obs + SSC * a[1]^2/2/a[2]
    ##else if (a[1] > 0)
    ##
    S.Chao1 <- S.obs + SSC * a[1]*(a[1]-1) / (a[2]+1)/2
    ##else
    ##    S.Chao1 <- S.obs
    Deriv.Ch1 <- gradF(a, i)

    ## The commonly used variance estimator is wrong for bias-reduced
    ## Chao estimate. It is based on the variance estimator of basic
    ## Chao estimate, but replaces the basic terms with corresponding
    ## terms in the bias-reduced estimate. The following is directly
    ## derived from the bias-reduced estimate.

    ## The commonly used one (for instance, in EstimateS):
    ##sd.Chao1 <-
    ##    sqrt(SSC*(a[1]*(a[1]-1)/2/(a[2]+1) +
    ##              SSC*(a[1]*(2*a[1]-1)^2/4/(a[2]+1)^2 +
    ##                  a[1]^2*a[2]*(a[1]-1)^2/4/(a[2]+1)^4)))

    sd.Chao1 <- (a[1]*((-a[2]^2+(-2*a[2]-a[1])*a[1])*a[1] +
                       (-1+(-4+(-5-2*a[2])*a[2])*a[2] +
                        (-2+(-1+(2*a[2]+2)*a[2])*a[2] +
                         (4+(6+4*a[2])*a[2] + a[1]*a[2])*a[1])*a[1])*S.Chao1))/
                             4/(a[2]+1)^4/S.Chao1
    sd.Chao1 <- sqrt(sd.Chao1)

    C.ace <- 1 - a[1]/N.rare
    i <- seq_along(a)
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

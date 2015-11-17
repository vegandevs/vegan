`mso` <-
    function (object.cca, object.xy, grain = 1, round.up = FALSE,
              permutations = 0) 
{
    EPS <- sqrt(.Machine$double.eps)
    if (inherits(object.cca, "mso")) {
        rm <- which(class(object.cca) == "mso")
        class(object.cca) <- class(object.cca)[-rm]
    }
    object <- object.cca
    xy <- object.xy
    N <- nrow(object$CA$Xbar)
    if (inherits(object, "rda")) 
        N <- 1
    ## we expect xy are coordinates and calculate distances, but a
    ## swift user may have supplied distances, and we use them.
    ## However, we won't test for distances in square matrices, but
    ## treat that as a user mistake and let it go.
    if (inherits(xy, "dist"))
        Dist <- xy
    else
        Dist <- dist(xy)
    object$grain <- grain
    if (round.up) 
        H <- ceiling(Dist/grain) * grain
    else H <- round(Dist/grain) * grain
    hmax <- round((max(Dist)/2)/grain) *grain
    H[H > hmax] <- max(H)
    object$H <- H
    H <- as.vector(H)
    Dist <- sapply(split(Dist, H), mean)
    object$vario <- data.frame(H = names(table(H)), Dist = Dist, 
                               n = as.numeric(table(H)))
    test <- list()
    if (is.numeric(object$CCA$rank)) {
        if (is.numeric(object$pCCA$rank)) {
            test$pcca <- sapply(split(dist(object$pCCA$Fit)^2 * 
                                      N/2, H), mean)
            test$cca <- sapply(split(dist(object$CCA$Xbar - object$CA$Xbar)^2 * 
                                     N/2, H), mean)
            test$ca <- sapply(split(dist(object$CA$Xbar)^2 * 
                                    N/2, H), mean)
            test$all.cond <- sapply(split(dist(object$CCA$Xbar)^2 * 
                                          N/2, H), mean)
            test$se <- sqrt(sapply(split(dist(object$CCA$Xbar)^2 * 
                                         N/2, H), var)/object$vario$n)
            object$vario <- cbind(object$vario, All = test$all.cond, 
                                  Sum = test$ca + test$cca, CA = test$ca,
                                  CCA = test$cca,  pCCA = test$pcca,
                                  se = test$se)
        } else {
            test$all <- sapply(split(dist(object$CCA$Xbar)^2 * 
                                     N/2, H), mean)
            test$cca <- sapply(split(dist(object$CCA$Xbar - object$CA$Xbar)^2 * 
                                     N/2, H), mean)
            test$ca <- sapply(split(dist(object$CA$Xbar)^2 * 
                                    N/2, H), mean)
            test$se <- sqrt(sapply(split(dist(object$CCA$Xbar)^2 * 
                                         N/2, H), var)/object$vario$n)
            object$vario <- cbind(object$vario, All = test$all, 
                                  Sum = test$ca + test$cca, CA = test$ca,
                                  CCA = test$cca, se = test$se)
        }
    } else {
        test$ca <- sapply(split(dist(object$CA$Xbar)^2 * N/2, 
                                H), mean)
        object$vario <- cbind(object$vario, All = test$ca, CA = test$ca)
    }
    permat <- getPermuteMatrix(permutations, nrow(object$CA$Xbar))
    nperm <- nrow(permat)
    if (nperm) {
        object$H.test <- matrix(0, length(object$H), nrow(object$vario))
        for (i in 1:nrow(object$vario)) {
            object$H.test[, i] <- as.numeric(object$H == object$vario$H[i])
        }
        xdis <- as.matrix(dist(object$CA$Xbar)^2)
        ## taking lower triangle is faster than as.dist() because it
        ## does not set attributes
        ltri <- lower.tri(xdis)
        statistic <- abs(cor(as.vector(xdis[ltri]), object$H.test))
        permfunc <- function(k) {
            permvec <- as.vector(xdis[k,k][ltri])
            abs(cor(permvec, object$H.test))
        }
        perm <- sapply(1:nperm, function(take) permfunc(permat[take,]))
        object$vario$CA.signif <-
            (rowSums(sweep(perm, 1, statistic - EPS, ">=")) + 1)/
                (nperm + 1)
        attr(object$vario, "control") <- attr(permat, "control")
    }
    object$call <- match.call()
    class(object) <- c("mso", class(object))
    object
}

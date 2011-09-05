`mso` <-
    function (object.cca, object.xy, grain = 1, round.up = FALSE,
              permutations = FALSE) 
{
    if (inherits(object.cca, "mso")) {
        rm <- which(class(object.cca) == "mso")
        class(object.cca) <- class(object.cca)[-rm]
    }
    object <- object.cca
    xy <- object.xy
    N <- nrow(object$CA$Xbar)
    if (inherits(object, "rda")) 
        N <- 1
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
    if (permutations) {
        ##require(base)
        object$H.test <- matrix(0, length(object$H), nrow(object$vario))
        for (i in 1:nrow(object$vario)) {
            object$H.test[, i] <- as.numeric(object$H == object$vario$H[i])
        }
        xdis <- dist(object$CA$Xbar)^2
        N <- attr(xdis, "Size")
        statistic <- abs(cor(as.vector(xdis), object$H.test))
        perm <- matrix(0, length(statistic), permutations)
        for (i in 1:permutations) {
            take <- sample(N, N)
            permvec <- as.vector(as.dist(as.matrix(xdis)[take, 
                                                         take]))
            perm[, i] <- abs(cor(permvec, object$H.test))
        }
        object$vario$CA.signif <- apply((perm >= matrix(statistic, 
                                         nrow(perm), ncol(perm)))/permutations, 1, sum)
    }
    object$call <- match.call()
    class(object) <- c("mso", class(object))
    object
}


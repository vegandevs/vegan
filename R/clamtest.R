## CLAM, reproduction of software described in Chazdon et al. 2011
## Ecology, 92, 1332--1343
clamtest <-
function(comm, groups, coverage.limit = 10,
specialization = 2/3, npoints = 20, alpha = 0.05/20)
{
    ## inital checks
    comm <- as.matrix(comm)
    if (NROW(comm) < 2)
        stop("'comm' must have at least two rows")
    if (nrow(comm) > 2 && missing(groups))
        stop("'groups' is missing")
    if (nrow(comm) == 2 && missing(groups))
        groups <- if (is.null(rownames(comm)))
            c("Group.1", "Group.2") else rownames(comm)
    if (length(groups) != nrow(comm))
        stop("length of 'groups' must equal to 'nrow(comm)'")
    if (length(unique(groups)) != 2)
        stop("number of groups must be two")
    glabel <- as.character(unique(groups))
    if (is.null(colnames(comm)))
        colnames(comm) <- paste("Species", 1:ncol(comm), sep=".")
    if (any(colSums(comm) <= 0))
        stop("'comm' contains columns with zero sums")
    spp <- colnames(comm)
    ## reproduced from Chazdon et al. 2011, Ecology 92, 1332--1343
    S <- ncol(comm)
    if (nrow(comm) == 2) {
        Y <- comm[glabel[1],]
        X <- comm[glabel[2],]
    } else {
        Y <- colSums(comm[which(groups==glabel[1]),])
        X <- colSums(comm[which(groups==glabel[2]),])
    }
    names(X) <- names(Y) <- NULL
    #all(ct$Total_SG == Y)
    #all(ct$Total_OG == X)
    m <- sum(Y)
    n <- sum(X)
    if (sum(Y) <= 0 || sum(X) <= 0)
        stop("group totals of zero are not allowed")
    ## check if comm contains integer, especially for singletons
    if (any(X[X>0] < 1) || any(Y[Y>0] < 1))
        warning("non-integer values <1 detected: analysis may not be meaningful")
    if (abs(sum(X,Y) - sum(as.integer(X), as.integer(Y))) > 10^-6)
        warning("non-integer values detected")
    C1 <- 1 - sum(X==1)/n
    C2 <- 1 - sum(Y==1)/m
    ## this stands for other than 2/3 cases
    uu <- specialization/(1-specialization)
    ## critical level
    Zp <- qnorm(alpha, lower.tail=FALSE)
    #p_i=a
    #pi_i=b
    ## function to calculate test statistic from Appendix D
    ## (Ecological Archives E092-112-A4)
    ## coverage limit is count, not freq !!!
    testfun <- function(p_i, pi_i, C1, C2, n, m) {
        C1 <- ifelse(p_i*n < coverage.limit, C1, 1)
        C2 <- ifelse(pi_i*m < coverage.limit, C2, 1)
        Var <- C1^2*(p_i*(1-p_i)/n) + uu^2*C2^2*(pi_i*(1-pi_i)/m)
        C1*p_i - C2*pi_i*uu - Zp*sqrt(Var)
    }
    ## root finding for iso-lines (instead of itarative search)
    rootfun <- function(pi_i, C1, C2, n, m, upper) {
        f <- function(p_i) testfun(p_i/n, pi_i/m, C1, C2, n, m)
        if (length(unique(sign(c(f(1), f(upper))))) > 1)
        ceiling(uniroot(f, lower=1, upper=upper)$root) else NA
    }
    ## sequences for finding Xmin and Ymin values
    Xseq <- as.integer(trunc(seq(1, max(X), len=npoints)))
    Yseq <- as.integer(trunc(seq(1, max(Y), len=npoints)))
    ## finding Xmin and Ymin values for Xseq and Yseq
    Xmins <- sapply(Yseq, function(z) rootfun(z, C1, C2, n, m, upper=max(X)))
    Ymins <- sapply(Xseq, function(z) rootfun(z, C2, C1, m, n, upper=max(Y)))

    ## needed to tweak original set of rules (extreme case reported
    ## by Richard Telford failed here)
    if (all(is.na(Xmins)))
        Xmins[1] <- 1
    if (all(is.na(Ymins)))
        Ymins[1] <- 1

    minval <- list(data.frame(x=Xseq[!is.na(Ymins)], y=Ymins[!is.na(Ymins)]),
        data.frame(x=Xmins[!is.na(Xmins)], y=Yseq[!is.na(Xmins)]))

    ## shared but too rare
    Ymin <- Ymins[1]
    Xmin <- Xmins[1]
    sr <- X < Xmin & Y < Ymin

    ## consequence of manually setting Xmin/Ymin resolved here
    tmp1 <- if (Xmin==1)
        list(x=1, y=Xmin) else approx(c(Xmin, 1), c(1, Ymin), xout=1:Xmin)
    tmp2 <- if (Ymin==1)
        list(x=1, y=Ymin) else approx(c(1, Ymin), c(Xmin, 1), xout=1:Ymin)

    for (i in 1:S) {
        if (X[i] %in% tmp1$x)
            sr[i] <- Y[i] < tmp1$y[which(X[i]==tmp1$x)]
        if (Y[i] %in% tmp2$x)
            sr[i] <- X[i] < tmp2$y[which(Y[i]==tmp2$x)]
    }
    ## classification
    a <- ifelse(X==0, 1, X)/n # \hat{p_i}
    b <- ifelse(Y==0, 1, Y)/m # \hat{\pi_i}
    specX <- !sr & testfun(a, b, C1, C2, n, m) > 0
    specY <- !sr & testfun(b, a, C2, C1, m, n) > 0
    gen <- !sr & !specX & !specY
    ## crosstable
    tmp <- ifelse(cbind(gen, specY, specX, sr), 1, 0)
    classes <- factor((1:4)[rowSums(tmp*col(tmp))], levels=1:4)
    levels(classes) <- c("Generalist", paste("Specialist", glabel[1], sep="_"),
        paste("Specialist", glabel[2], sep="_"), "Too_rare")
    tab <- data.frame(Species=spp, y=Y, x=X, Classes=classes)
    colnames(tab)[2:3] <- paste("Total", glabel, sep="_")
    rownames(tab) <- NULL
    class(tab) <- c("clamtest","data.frame")
    attr(tab, "settings") <- list(labels = glabel,
        coverage.limit = coverage.limit, specialization = specialization,
        npoints = npoints, alpha = alpha)
    attr(tab, "minv") <- minval
    attr(tab, "coverage") <- structure(c(C2, C1), .Names=glabel)
    tab
}

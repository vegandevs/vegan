## adipart function START
adipart <-
function(matr, strata, hclass=NULL, method="trad", index=c("richness", "shannon", "simpson"), 
    scales=seq(0, 2, 0.2), weights="unif", test=TRUE, permtype="full", times=100, crit=0.05, 
    burnin=10000, results=FALSE, ...)
{
    method <- match.arg(method, c("trad", "tsallis"))
    if (method == "trad")
        index <- match.arg(index, c("richness", "shannon", "simpson"), is.null(hclass))
    weights <- match.arg(weights, c("unif", "prop"))
    permtype <- match.arg(permtype, c("full", "swap"))
    if (length(unique(strata[,1])) != 1)
        stop("first column of strata should be uniform")
    if (method == "tsallis" && length(scales) != 1 && !is.null(hclass))
        stop("scales must be of length 1 if hclass is defined")
    if (times == 0) test <- FALSE
    if (!test && results)
        stop("results are not produced when test is FALSE or times = 0")
    if (!is.null(hclass) && results)
        stop("results are not produced when hclass is defined")
    if (inherits(matr, "matrix") || inherits(matr, "data.frame")) m <- matr
    if (inherits(matr, "permat")) m <- matr$orig
    if (test) {
        if (inherits(matr, "permat")) {
            perm <- matr$perm
            times <- attr(matr, "times")
            perm.done <- TRUE
            } else {
                 perm <- NULL
                 perm.done <- FALSE}
        } else perm <- NULL
## internal
nestFactor <- function(stra) {
    nr <- nrow(stra)
    nc <- ncol(stra)
    for (k in 1:nc) stra[,k] <- as.numeric(stra[,k])
    out <- stra
    stra <- data.frame(rep(1,nrow(stra)), stra)
    for (i in 1:nc) {
        if (i == 1) out[,i] <- interaction(stra[,i], stra[,(i+1)], drop=TRUE)
        else out[,i] <- interaction(out[,(i-1)], stra[,(i+1)], drop=TRUE)
        levels(out[,i]) <- 1:nlevels(out[,i])}
    out <- out[,c(nc:1)]
    colnames(out) <- paste("x", 1:ncol(out), sep="")
    rownames(out) <- 1:nrow(out)
    return(out)}
## internal !!!f <- nestFactor(x)
adpTrad <- function(y, f, index, weights="unif", serr=TRUE){
    if (any(rowSums(y) == 0)) stop("empty samples not allowed")
    ni <- length(index)
    n <- ncol(f)
    q <- list()
    w <- list()
    a <- list()
    le <- list()
    if (serr) sde <- list()
    for (i in 1:n){
        tab <- aggregate(y, by=list(f[,i]), sum)[,-1]
        rndq <- character(3)
        if ("richness" %in% index) {
            Richness <- apply(tab > 0, 1, sum)
            rndq[1] <- "Richness"
            } else {Richness <- NULL
            rndq[1] <- "X"}
        if ("shannon" %in% index) {
            Shannon <- diversity(tab, "shannon")
            rndq[2] <- "Shannon"
            } else {Shannon <- NULL
            rndq[2] <- "X"}
        if ("simpson" %in% index) {
            Simpson <- diversity(tab, "simpson")
            rndq[3] <- "Simpson"
            } else {Simpson <- NULL
            rndq[3] <- "X"}
        rndq <- rndq[which(rndq != "X")]
        dq.list <- list(Richness, Shannon, Simpson)
        dq <- dqsd <- data.frame(matrix(unlist(dq.list), nrow(tab), length(rndq)))
        colnames(dq) <- rndq 
#        dq[dq == Inf | dq == -Inf] <- 0
        le[[i]] <- length(unique(f[,i]))
        if (weights=="prop")
            q[[i]] <- tapply(apply(y, 1, sum), list(f[,i]), sum) / sum(y)
        if (weights=="unif")
            q[[i]] <- rep(1 / length(unique(f[,i])), length(unique(f[,i])))
        if (ncol(as.matrix(dq)) != 1) {
            for (j in 1:ni) dq[,j] <- as.matrix(dq)[,j] * q[[i]]
            } else {dq <- dq * q[[i]]}
        w[[i]] <- dq
        if (ni == 1) {
            a[[i]] <- sum(w[[i]])
            if (serr) sde[[i]] <- sd(dqsd)
        } else {
        if (ncol(as.matrix(dq)) != 1) {
            a[[i]] <- apply(as.matrix(w[[i]]), 2, sum)
            if (serr) sde[[i]] <- apply(as.matrix(dqsd), 2, sd)
            } else {
            a[[i]] <- w[[i]]
            if (serr) sde[[i]] <- sd(dqsd)}}
        }
    alpha <- matrix(unlist(a), n, ni, byrow=TRUE)
    beta <- matrix(alpha[1:(n-1),], (n-1), ni, byrow=TRUE)
    for (i in 1:(n-1))
        beta[i, 1:ni] <- alpha[(i+1), 1:ni] - alpha[i, 1:ni]
    colnames(alpha) <- colnames(beta) <- rndq
    if (serr) {
        sdev <- matrix(unlist(sde), n, ni, byrow=TRUE)
        sqrn <- matrix(unlist(le), n, ni, byrow=FALSE)
        sem <- sdev / sqrt(sqrn)
#        sem[is.na(sem)] <- 0
        colnames(sem) <- rndq}
    if (serr) {out <- list(alpha=alpha, beta=beta, sem=sem)
        } else out <- list(alpha=alpha, beta=beta)
    return(out)}
## internal !!!f <- nestFactor(x)
adpTsallis <- function(y, f, weights="unif", serr=TRUE, scales=seq(0, 2, 0.2)){
    if (any(rowSums(y) == 0))
        stop("empty samples not allowed")
    n <- ncol(f)
    ni <- length(scales)
    q <- list()
    w <- list()
    a <- list()
    le <- list()
    if (serr) sde <- list()
    for (i in 1:n){
        tab <- aggregate(y, by=list(f[,i]), sum)[,-1]
        dq <- tsallis(tab, scales=scales, norm=FALSE)
        if (i == 1) rndq <- colnames(dq)
        if (ni == 1) dq <- data.frame(matrix(dq, length(dq), 1))
        dqsd <- dq
#        dq[dq == Inf | dq == -Inf] <- 0
        le[[i]] <- length(unique(f[,i]))
        if (weights=="prop")
            q[[i]] <- tapply(apply(y, 1, sum), list(f[,i]), sum) / sum(y)
        if (weights=="unif")
            q[[i]] <- rep(1 / length(unique(f[,i])), length(unique(f[,i])))
        if (ncol(as.matrix(dq)) != 1) {
            for (j in 1:ni) dq[,j] <- as.matrix(dq)[,j] * q[[i]]
            } else {dq <- dq * q[[i]]}
        w[[i]] <- dq
        if (any(dim(w[[i]]) != 1)) {
            a[[i]] <- apply(w[[i]], 2, sum)
            if (serr) sde[[i]] <- apply(dqsd, 2, sd)
            } else {
            if (ni == 1) {
                a[[i]] <- apply(w[[i]], 1, sum)
                if (serr) sde[[i]] <- sd(dqsd)}
            if (i == n) {
                a[[i]] <- w[[i]]
                if (serr) sde[[i]] <- apply(as.matrix(dqsd), 1, sd)}
            }
        }
    alpha <- matrix(unlist(a), n, ni, byrow=TRUE)
    beta <- matrix(alpha[1:(n-1),], (n-1), ni, byrow=TRUE)
    for (i in 1:(n-1))
        beta[i, 1:ni] <- alpha[(i+1), 1:ni] - alpha[i, 1:ni]
    colnames(alpha) <- colnames(beta) <- rndq
    if (serr) {
        sdev <- matrix(unlist(sde), n, ni, byrow=TRUE)
        sqrn <- matrix(unlist(le), n, ni, byrow=FALSE)
        sem <- sdev / sqrt(sqrn)
#        sem[is.na(sem)] <- 0
        colnames(sem) <- rndq}
    if (serr) {out <- list(alpha=alpha, beta=beta, sem=sem)
        } else out <- list(alpha=alpha, beta=beta)
    return(out)}
## internal
formatAdp <- function(x, style=c("alpha","beta","alpha2","beta2", "alphaH"),col.nam=NULL){
    if (style=="alpha" || style=="alpha2" || style=="alphaH") {
        rownames(x) <- interaction(rep("Alpha",nrow(x)), 1:nrow(x))
        if (style=="alpha" || style=="alphaH") rownames(x)[nrow(x)] <- "Gamma"
        if (style=="alpha2" || style=="alphaH") colnames(x) <- col.nam}
    if (style=="beta" || style=="beta2") {
        rownames(x) <- interaction(rep("Beta",nrow(x)), 1:nrow(x))
        if (style=="beta2") colnames(x) <- col.nam}
    return(x)}
## internal
adpPvalue <- function(obs, perm, crit=0.05, style="alpha2", col.nam=NULL){
#    perm[perm < 0] <- 0
    Mean <- apply(perm, 1, mean)
    obs2 <- Pval <- array(obs)
    Sign <- (Mean - obs2) > 0
    Abs <- abs(Mean - obs2)
    sdp <- apply(perm, 1, sd)
    upper <- lower <- p.val <- cl1 <- cl2 <- numeric(length(obs2))
    adp1 <- Mean + Abs
    adp2 <- Mean - Abs
    for (i in 1:(nrow(obs)*ncol(obs))) {
        upper[i] <- sum(perm[i,] >= adp1[i])
        lower[i] <- sum(perm[i,] <= adp2[i])
        p.val[i] <- if (Sign[i]) upper[i] else lower[i]
        cl1[i] <- quantile(perm[i,], probs=crit/2)
        cl2[i] <- quantile(perm[i,], probs=1-crit/2)}
    Pval <- matrix(p.val*2 / times, nrow(obs), ncol(obs))
#    Pval[Pval > 1] <- 1
    Mean <- matrix(Mean, nrow(obs), ncol(obs))
    cl1 <- matrix(cl1, nrow(obs), ncol(obs))
    cl2 <- matrix(cl2, nrow(obs), ncol(obs))
    ses <- (obs - Mean) / matrix(sdp, nrow(obs), ncol(obs))
    Pval <- formatAdp(Pval, style, col.nam)
    Mean <- formatAdp(Mean, style, col.nam)
    cl1 <- formatAdp(cl1, style, col.nam)
    cl2 <- formatAdp(cl2, style, col.nam)
    ses <- formatAdp(ses, style, col.nam)
    return(list(p.value=Pval,mean=Mean,cl1=cl1,cl2=cl2,ses=ses))}
## internal fun
simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")}
## internal
testAdp <- function(obs, perm, f, method, results, type=c(1,2), id=NULL, burnin=NULL, ...){
#    if (type==1) id <- 1:nrow(perm[[1]])
    pa <- matrix(NA, (nrow(obs$alpha)-1)*ncol(obs$alpha), times)
    pb <- matrix(NA, nrow(obs$beta)*ncol(obs$alpha), times)
    for (i in 1:times) {
        if (perm.done) {
            perm.i <- perm[[i]]
            } else {
                if (permtype == "full")
                    perm.i <- permatfull(matr, times=1, ...)$perm[[1]]
                if (permtype == "swap") {
                    if (i == 1) {
                        perm.i <- permatswap(matr, times=1, burnin=burnin, method="swap", ...)$perm[[1]]
                        } else {
                        perm.i <- permatswap(perm.i, times=1, burnin=0, method="swap", ...)$perm[[1]]
                }}}
        if (method == "trad")
            adp.perm <- adpTrad(perm.i[id,], f[id,], index, weights, serr=FALSE)
        if (method == "tsallis")
            adp.perm <- adpTsallis(perm.i[id,], f[id,], weights, serr=FALSE, scales=scales)
        pa[,i] <- array(adp.perm$alpha[1:(nrow(adp.perm$alpha)-1),])
        pb[,i] <- array(adp.perm$beta)}
    obs.in.a <- obs$alpha[1:(nrow(adp.perm$alpha)-1),]
    if (any(dim(as.matrix(obs.in.a)) == 1))
        obs.in.a <- matrix(obs.in.a, length(obs.in.a),1)
    obs.in.b <- obs$beta
    if (any(dim(as.matrix(obs.in.b)) == 1))
    obs.in.a <- matrix(obs.in.b, length(obs.in.b),1)
    test.alpha <- adpPvalue(obs.in.a,pa,crit,"alpha2",colnames(obs$alpha))
    test.beta <- adpPvalue(obs.in.b,pb,crit,"beta2",colnames(obs$beta))
    if (results) res <- list(alpha=t(pa), beta=t(pb)) else res <- NULL
    out <- list(alpha=test.alpha, beta=test.beta, res=res)
    return(out)}

## make orig
    if (!is.null(hclass)) {
        habnam <- unique(as.character(hclass))
        nhab <- length(habnam)
        habnam2 <- c(habnam, "All")
        f2 <- data.frame(strata[,1], hclass, strata[,2:ncol(strata)])
        obslist <- list(alpha=list(), beta=list(), sem=list())
        for (i in 1:(nhab+1)) {
            if (i < nhab+1) {
                ff <- nestFactor(f2[hclass == habnam[i],])
                mm <- m[hclass == habnam[i],]
            } else {
                ff <- nestFactor(f2)
                mm <- m}
            if (method == "trad")
                obs <- adpTrad(mm, ff, index, weights, TRUE)
            if (method == "tsallis")
                obs <- adpTsallis(mm, ff, weights, TRUE, scales=scales)
            obslist$alpha[[i]] <- obs$alpha
            obslist$beta[[i]] <- obs$beta
            obslist$sem[[i]] <- obs$sem
        }
        alpha <- formatAdp(matrix(unlist(obslist$alpha), ncol(f2), nhab+1),"alphaH",habnam2)
        beta <- formatAdp(matrix(unlist(obslist$beta), ncol(f2)-1, nhab+1),"beta2",habnam2)
        sem <- formatAdp(matrix(unlist(obslist$sem), ncol(f2), nhab+1),"alphaH",habnam2)
        obs <- list(alpha=alpha, beta=beta, sem=sem)
        if (test) {
            explist <- list(
                alpha=list(p.value=list(), mean=list(), cl1=list(), cl2=list(), ses=list()),
                beta=list(p.value=list(), mean=list(), cl1=list(), cl2=list(), ses=list()))
            for (i in 1:(nhab+1)) {
                if (i < nhab+1) {
                    expt <- testAdp(obs, perm, ff, method, FALSE, 2, which(hclass == habnam[i]), burnin, ...)
                } else {
                    expt <- testAdp(obs, perm, ff, method, FALSE, 2, 1:length(hclass), burnin, ...)
                }
            explist$alpha$p.value[[i]] <- expt$alpha$p.value
            explist$alpha$mean[[i]] <- expt$alpha$mean
            explist$alpha$cl1[[i]] <- expt$alpha$cl1
            explist$alpha$cl2[[i]] <- expt$alpha$cl2
            explist$alpha$ses[[i]] <- expt$alpha$ses
            explist$beta$p.value[[i]] <- expt$beta$p.value
            explist$beta$mean[[i]] <- expt$beta$mean
            explist$beta$cl1[[i]] <- expt$beta$cl1
            explist$beta$cl2[[i]] <- expt$beta$cl2
            explist$beta$ses[[i]] <- expt$beta$ses
        }
        alpha.p.value <- formatAdp(matrix(unlist(explist$alpha$p.value), ncol(f2)-1, nhab+1),"alpha2",habnam2)
        alpha.mean <- formatAdp(matrix(unlist(explist$alpha$mean), ncol(f2)-1, nhab+1),"alpha2",habnam2)
        alpha.cl1 <- formatAdp(matrix(unlist(explist$alpha$cl1), ncol(f2)-1, nhab+1),"alpha2",habnam2)
        alpha.cl2 <- formatAdp(matrix(unlist(explist$alpha$cl2), ncol(f2)-1, nhab+1),"alpha2",habnam2)
        alpha.ses <- formatAdp(matrix(unlist(explist$alpha$ses), ncol(f2)-1, nhab+1),"alpha2",habnam2)
        beta.p.value <- formatAdp(matrix(unlist(explist$beta$p.value), ncol(f2)-1, nhab+1),"beta2",habnam2)
        beta.mean <- formatAdp(matrix(unlist(explist$beta$mean), ncol(f2)-1, nhab+1),"beta2",habnam2)
        beta.cl1 <- formatAdp(matrix(unlist(explist$beta$cl1), ncol(f2)-1, nhab+1),"beta2",habnam2)
        beta.cl2 <- formatAdp(matrix(unlist(explist$beta$cl2), ncol(f2)-1, nhab+1),"beta2",habnam2)
        beta.ses <- formatAdp(matrix(unlist(explist$beta$ses), ncol(f2)-1, nhab+1),"beta2",habnam2)
        exp <- list(
            alpha=list(p.value=alpha.p.value, mean=alpha.mean, cl1=alpha.cl1, cl2=alpha.cl2, ses=alpha.ses),
            beta=list(p.value=beta.p.value, mean=beta.mean, cl1=beta.cl1, cl2=beta.cl2, ses=beta.ses))
        res <- NULL
        } else {
            exp <- NULL
            times <- 0
            res <- NULL}
        obs$alpha[(nrow(obs$alpha)-1), 1:(ncol(obs$alpha)-1)] <- NA
        obs$beta[nrow(obs$beta), 1:(ncol(obs$beta)-1)] <- NA
        if (test) {for (a in 1:5) {
            exp$alpha[[a]][(nrow(obs$alpha)-1), 1:(ncol(obs$alpha)-1)] <- NA
            exp$beta[[a]][(nrow(obs$alpha)-1), 1:(ncol(obs$alpha)-1)] <- NA
            }}
        index.out <- simpleCap(index)
        design <- "twoway"
    } else {
        f <- nestFactor(strata)
        if (method == "trad")
            obs <- adpTrad(m, f, index, weights, TRUE)
        if (method == "tsallis")
            obs <- adpTsallis(m, f, weights, TRUE, scales=scales)
        obs$alpha <- formatAdp(obs$alpha,"alpha")
        obs$sem <- formatAdp(obs$sem,"alpha")
        obs$beta <- formatAdp(obs$beta,"beta")
        index.out <- colnames(obs$alpha)
        design <- "oneway"
        if (test) {
            exp <- testAdp(obs, perm, f, method, results, 1, 1:nrow(m), burnin, ...)
            res <- exp$res
            exp <- exp[1:2]
        } else {
            exp <- NULL
            times <- 0
            res <- NULL}
    }

    input <- list(call=match.call(), m=m, f=nestFactor(strata), h=hclass)
    out <- list(input=input, obs=obs, exp=exp, res=res)
    attr(out, "method") <- method
    attr(out, "design") <- design
    if (method == "trad") attr(out, "index") <- index.out
    if (method == "tsallis") attr(out, "index") <- scales
    attr(out, "times") <- times
    attr(out, "weights") <- weights
    attr(out, "crit") <- crit
    class(out) <- c("adipart", "list")
    return(out)
} ## adp function END

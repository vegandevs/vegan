eventstar <- function(x, qmax=5) {
    if (is.null(dim(x)))
        x <- matrix(x, 1, length(x))
    lossfun <- function(q, x)
        tsallis(x, scales=q, norm=TRUE)
    qstarfun <- function(x) {
        optimize(lossfun, interval=c(0, qmax), x=x)$minimum
    }
    qs <- apply(x, 1, qstarfun)
    Hs <- sapply(1:nrow(x), function(i) tsallis(x[i,], 
        scales=qs[i], hill=FALSE))
    S <- rowSums(x)
    Es <- ifelse(qs==1, log(S), Hs/((S^(1-qs)-1)/(1-qs)))
    Ds <- (1-(qs-1)*Hs)^(1/(1-qs))
    out <- data.frame(qstar=qs, Estar=Es, Hstar=Hs, Dstar=Ds)
    rownames(out) <- rownames(x)
    out
}

`taxondive` <-
    function (comm, dis, match.force = FALSE)
{
    binary <- FALSE
    comm <- as.matrix(comm)
    if (missing(dis)) {
        n <- ncol(comm)
        dis <- structure(rep(1, n * (n - 1)/2), Size = n, class = "dist")
    }
    dis <- as.dist(dis)
    if (match.force || attr(dis, "Size") != ncol(comm)) {
        if (match.force)
            message("forced matching 'dis' labels and 'comm' names")
        else
            message("dimensions do not match between 'comm' and 'dis'")
        if (all(colnames(comm) %in% labels(dis))) {
            dis <- as.matrix(dis)
            dis <- as.dist(dis[colnames(comm), colnames(comm)])
            message("matched 'dis' labels by 'comm' names")
        } else {
            stop("could not match names in 'dis' and 'comm'")
        }
        if (length(unique(colnames(comm))) != ncol(comm))
            stop("names not in unique in 'comm': match wrong")
        if (length(unique(labels(dis))) != attr(dis, "Size"))
            warning("labels not unique in 'dis': matching probably wrong")
    }
    del <- dstar <- dplus <- Ed <- Edstar <- edplus <- NULL
    if (!binary) {
        del <- apply(comm, 1, function(x) sum(as.dist(outer(x,
                                                            x)) * dis))
        dstar <- apply(comm, 1, function(x) sum(dis * (xx <- as.dist(outer(x,
                                                                           x))))/sum(xx))
        rs <- rowSums(comm)
        del <- del/rs/(rs - 1) * 2
        cs <- colSums(comm)
        tmp <- sum(as.dist(outer(cs, cs)) * dis)
        Ed <- tmp/sum(cs)/sum(cs - 1) * 2
        Edstar <- tmp/sum(cs)/(sum(cs) - 1) * 2
    }
    comm <- ifelse(comm > 0, 1, 0)
    dplus <- apply(comm, 1, function(x) sum(as.dist(outer(x,
                                                          x)) * dis))
    Lambda <- apply(comm, 1, function(x) sum(as.dist(outer(x,
                                                           x)) * dis^2))
    m <- rowSums(comm)
    dplus <- dplus/m/(m - 1) * 2
    Lambda <- Lambda/m/(m - 1) * 2 - dplus^2
    S <- attr(dis, "Size")
    omebar <- sum(dis)/S/(S - 1) * 2
    varome <- sum(dis^2)/S/(S - 1) * 2 - omebar^2
    omei <- rowSums(as.matrix(dis))/(S - 1)
    varomebar <- sum(omei^2)/S - omebar^2
    vardplus <- 2 * (S - m)/(m * (m - 1) * (S - 2) * (S - 3)) *
        ((S - m - 1) * varome + 2 * (S - 1) * (m - 2) * varomebar)
    out <- list(Species = m, D = del, Dstar = dstar, Lambda = Lambda,
                Dplus = dplus, sd.Dplus = sqrt(vardplus), SDplus = m *
                dplus, ED = Ed, EDstar = Edstar, EDplus = omebar)
    class(out) <- "taxondive"
    out
}

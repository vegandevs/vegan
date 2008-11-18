"nestedn0" <-
function(comm)
{
    comm <- ifelse(comm > 0, 1, 0)
    R <- rowSums(comm)
    spmin <- apply(comm, 2, function(x) min((x*R)[x > 0]))
    n0 <- spmin
    for (i in 1:ncol(comm))
        n0[i] <- sum(comm[,i] == 0 & R > spmin[i])
    out <- list(spmin = spmin, n0 = n0, statistic = sum(n0))
    class(out) <- "nestedn0"
    out
}


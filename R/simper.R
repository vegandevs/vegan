`simper` <-
    function(comm, group)
{
    comp <- t(combn(unique(as.character(group)), 2))
    outlist <- NULL
    for (i in 1:nrow(comp)){
        group.a <- as.matrix(comm[group == comp[i, 1], ])  
        group.b <- as.matrix(comm[group == comp[i, 2], ])  
        n.a <- nrow(group.a)
        n.b <- nrow(group.b)
        P <- ncol(comm)
        contr <- matrix(ncol = P, nrow = n.a * n.b)
        for(j in 1:n.b) {
            for(k in 1:n.a) {
                md <- abs(group.a[k, ] - group.b[j, ])
                me <- group.a[k, ] + group.b[j, ]
                contr[(j-1)*n.a+k, ] <- md / sum(me)	
            }
        }
        av.contr <- colMeans(contr) * 100
        ov.av.dis <- sum(av.contr)
        sdi <- apply(contr, 2, sd)
        sdi.av <- av.contr / sdi
        av.a <- colMeans(group.a)
        av.b <- colMeans(group.b) 
        ord <- order(av.contr, decreasing = TRUE)
        out <-  list(species = colnames(comm), average = av.contr, overall = ov.av.dis, sd = sdi, meansdratio = sdi.av, ava = av.a, avb = av.b, ord = ord)
        outlist[[paste(comp[i,1], "_", comp[i,2], sep = "")]] <- out
    }
    class(outlist) <- "simper"
    outlist
}

`print.simper` <-
    function(x, ...)
{
    cusum <- lapply(x, function(z) cumsum(z$average[z$ord] / z$overall * 100))
    spec <- lapply(x, function(z) z$species[z$ord])
    for(i in 1:length(cusum)) {
        names(cusum[[i]]) <- spec[[i]]
    }
    out <- lapply(cusum, function(z) z[z <= 70])
    print(out)
    invisible(x)
}

`summary.simper` <-
    function(object, ...)
{
    cusum <- lapply(object, function(z) cumsum(z$average[z$ord] / z$overall * 100))
    out <- lapply(object, function(z) data.frame(contr = z$average, sd = z$sd, 'contr/sd' = z$meansdratio, av.a = z$ava, av.b = z$avb)[z$ord, ])
    for(i in 1:length(out)) {
        out[[i]]$cum <- cusum[[i]]
    }
    class(out) <- "summary.simper"
    out
}

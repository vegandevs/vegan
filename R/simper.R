`simper` <-
    function(comm, group, permutations = 0)
{
    comm <- as.matrix(comm)
    comp <- t(combn(unique(as.character(group)), 2))
    outlist <- NULL
    for (i in 1:nrow(comp)) {
        group.a <- comm[group == comp[i, 1], ]
        group.b <- comm[group == comp[i, 2], ]
        n.a <- nrow(group.a)
        n.b <- nrow(group.b)
        P <- ncol(comm)
        contr <- matrix(ncol = P, nrow = n.a * n.b)
        for (j in 1:n.b) {
            for (k in 1:n.a) {
                md <- abs(group.a[k, ] - group.b[j, ])
                me <- group.a[k, ] + group.b[j, ]
                contr[(j-1)*n.a+k, ] <- md / sum(me)	
            }
        }
        average <- colMeans(contr) * 100
        
        if(permutations != 0){
            cat("Permuting", paste(comp[i,1], "_", comp[i,2], sep = ""), "\n")
            nobs <- length(group)
            perm.contr <- matrix(nrow=P, ncol=permutations)
            contrp <- matrix(ncol = P, nrow = n.a * n.b)
            for(p in 1:permutations){
                perm <- shuffle(nobs)
                groupp <- group[perm]
                ga <- comm[groupp == comp[i, 1], ] 
                gb <- comm[groupp == comp[i, 2], ]
                for(j in 1:n.b) {
                    for(k in 1:n.a) {
                        mdp <- abs(ga[k, ] - gb[j, ])
                        mep <- ga[k, ] + gb[j, ]
                        contrp[(j-1)*n.a+k, ] <- mdp / sum(mep)  
                    }
                }
                perm.contr[ ,p] <- colMeans(contrp) * 100
            }
        p <- (apply(apply(perm.contr, 2, function(x) x >= average), 1, sum) + 1) / (permutations + 1)
        } 
        else {
          p <- NULL
        }
        
        overall <- sum(average)
        sdi <- apply(contr, 2, sd)
        ratio <- average / sdi
        av.a <- colMeans(group.a)
        av.b <- colMeans(group.b) 
        ord <- order(average, decreasing = TRUE)
        cusum <- cumsum(average[ord] / overall * 100)
        out <-  list(species = colnames(comm), average = average, overall = overall, sd = sdi, ratio = ratio, ava = av.a, avb = av.b, ord = ord, cusum = cusum, p = p)
        outlist[[paste(comp[i,1], "_", comp[i,2], sep = "")]] <- out
    }
    class(outlist) <- "simper"
    outlist
}

`print.simper` <-
    function(x, ...)
{
    cat("cumulative contributions of most influential species:\n\n")
    cusum <- lapply(x, function(z) z$cusum)
    spec <- lapply(x, function(z) z$species[z$ord])
    for (i in 1:length(cusum)) {
        names(cusum[[i]]) <- spec[[i]]
    }
    out <- lapply(cusum, function(z) z[z <= 70])
    print(out)
    invisible(x)
}

`summary.simper` <-
    function(object, ordered = TRUE, ...)
{
    if (ordered == TRUE) {
        out <- lapply(object, function(z) data.frame(contr = z$average, sd = z$sd, ratio = z$ratio, av.a = z$ava, av.b = z$avb)[z$ord, ])
        cusum <- lapply(object, function(z) z$cusum)
        for(i in 1:length(out)) {
            out[[i]]$cumsum <- cusum[[i]]
            if(!is.null(object[[i]]$p)) {
                out[[i]]$p <- object[[i]]$p[object[[i]]$ord]
            }
        } 
    } 
    else {
        out <- lapply(object, function(z) data.frame(contr = z$average, sd = z$sd, 'contr/sd' = z$ratio, av.a = z$ava, av.b = z$avb, p = z$p))
    }
    class(out) <- "summary.simper"
    out
}

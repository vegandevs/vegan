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
        dat <- data.frame(av.contr, sdi, sdi.av, av.a, av.b)
        dat <- dat[order(dat$av.contr, decreasing = TRUE), ]
        cum <-  cumsum(dat$av.contr / ov.av.dis) * 100
        out <-  data.frame(dat, cum)
        names(out) <- c("contr", "sd", "contr/sd", paste("av_", comp[i, 1], sep = ""), paste("av_", comp[i, 2], sep = ""), "cum")
        outlist[[paste(comp[i,1], "_", comp[i,2], sep = "")]] <- out
    }
    class(outlist) <- "simper"
    outlist
}

`print.simper` <-
    function(x, ...)
{
    out <- lapply(x, function(z) t(z[z$cum <= 70 ,"cum", drop = FALSE]))
    print(out)
    invisible(x)
}

`summary.simper` <-
    function(object, ...)
{
    class(object) <- "summary.simper"
    object
}

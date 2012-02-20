`simper` <-
    function(comm, group, permutations = 0, trace = FALSE,  ...)
{
    comm <- as.matrix(comm)
    comp <- t(combn(unique(as.character(group)), 2))
    outlist <- NULL
    ## data parameters
    P <- ncol(comm)
    nobs <- nrow(comm)
    ## Make permutation matrix
    if (length(permutations) == 1) {
        perm <- shuffleSet(nobs, permutations, ...)
    } else {  # permutations is a matrix
        perm <- permutations
    }
    ## check dims (especially if permutations was a matrix)
    if (ncol(perm) != nobs)
        stop(gettextf("'permutations' have %d columns, but data have %d rows",
                          ncol(perm), nobs))
    ## OK: take number of permutations
    nperm <- nrow(perm)
    if (nperm > 0)
        perm.contr <- matrix(nrow=P, ncol=nperm)
    for (i in 1:nrow(comp)) {
        group.a <- comm[group == comp[i, 1], ]
        group.b <- comm[group == comp[i, 2], ]
        n.a <- nrow(group.a)
        n.b <- nrow(group.b)
        contr <- matrix(ncol = P, nrow = n.a * n.b)
        for (j in 1:n.b) {
            for (k in 1:n.a) {
                md <- abs(group.a[k, ] - group.b[j, ])
                me <- group.a[k, ] + group.b[j, ]
                contr[(j-1)*n.a+k, ] <- md / sum(me)	
            }
        }
        average <- colMeans(contr) * 100
        
        if(nperm > 0){
            if (trace)
                cat("Permuting", paste(comp[i,1], comp[i,2], sep = "_"), "\n")
            contrp <- matrix(ncol = P, nrow = n.a * n.b)
            for(p in 1:nperm){
                groupp <- group[perm[p,]]
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
        out <- list(species = colnames(comm), average = average,
                    overall = overall, sd = sdi, ratio = ratio, ava = av.a,
                    avb = av.b, ord = ord, cusum = cusum, p = p)
        outlist[[paste(comp[i,1], "_", comp[i,2], sep = "")]] <- out
    }
    attr(outlist, "permutations") <- nperm
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
    ## this probably fails with empty or identical groups that have 0/0 = NaN
    out <- lapply(cusum, function(z) z[seq_len(min(which(z >= 70)))])
    print(out)
    invisible(x)
}

`summary.simper` <-
    function(object, ordered = TRUE, digits = max(3, getOption("digits") - 3), ...)
{
    if (ordered) {
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
        out <- lapply(object, function(z) data.frame(cbind(contr = z$average, sd = z$sd, 'contr/sd' = z$ratio, av.a = z$ava, av.b = z$avb, p = z$p)))
    }
    attr(out, "digits") <- digits
    attr(out, "permutations") <- attr(object, "permutations")
    class(out) <- "summary.simper"
    out
}

`print.summary.simper`<-
    function(x, digits = attr(x, "digits"), ...)
{
    signif.stars <- getOption("show.signif.stars") && attr(x, "permutations") > 0
    starprint <- function(z) {
        if (signif.stars && any(z$p < 0.1)) {
            stars <- symnum(z$p, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                            symbols = c("***", "**", "*", ".", " "))
            z <- cbind(z, " " = format(stars))
        }
        z
    }
    out <- lapply(x, starprint)
    for (nm in names(out)) {
        cat("\nContrast:", nm, "\n\n")
        print(out[[nm]], digits = digits, ...)
    }
    if (signif.stars && any(sapply(x, function(z) z$p) < 0.1)) {
        leg <- attr(symnum(1, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                            symbols = c("***", "**", "*", ".", " ")), "legend")
        cat("---\nSignif. codes: ", leg, "\n")
    }
    if ((np <- attr(x, "permutations")) > 0)
        cat("P-values based on", np, "permutations\n")
    invisible(x)
}

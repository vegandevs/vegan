#' Similarity Percentages
#'
#' Discriminating species between two groups using Bray-Curtis dissimilarities.
#'
#' @param comm Community data matrix
#' @param group factor describing the group structure. Must have at least 2 levels.
#' 
#' @section Details: simper is a function that extracts species discriminating between two groups and ranks them by their contribution. They Bray-Curtis dissimilarity between two samples is given as:
#' d[jk] = (sum abs(x[ij]-x[ik])) / (sum (x[ij]+x[ik]))
#' 
#' print.simper prints the species with a cumulative contribution <= 70\%
#' 
#' summary.simper shows the whole output with additional characteristics
#'
#' @return
#' A list of dataframes for every factor-combination.
#' \item{contr}{average contribution to overall dissimilarity}
#' \item{sd}{standard deviation of contribution}
#' \item{contr/sd}{mean to sd ratio}
#' \item{av_}{average abundance per group} 
#' \item{cum}{cumulative per cent contribution}
#'
#' @author Eduard Szöcs \email{szoe8822@@uni-landau.de}
#' @references Clarke, K.R. 1993. Non-parametric multivariate analyses of changes incommunity structure. Austral Ecology, 18(1):117–143.
#' @export
#'
#' @keywords multivariate
#'
#' @examples
#' require(vegan)
#' data(dune)
#' data(dune.env)
#' (sim <- simper(dune, dune.env$Management))
#' summary(sim)

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
        me <- matrix(ncol = P)    	
        md <- matrix(ncol = P)
        contr <- matrix(ncol = P, nrow = n.a * n.b)
        for(j in 1:n.b) {
            for(k in 1:n.a) {
                for(s in 1:P) {
                    md[s] <- abs(group.a[k, s] - group.b[j, s])
                    me[s] <- group.a[k, s] + group.b[j, s]
                    contr[(j-1)*n.a+k, ] <- md / rowSums(me)	
                }
            }
        }
        av.contr <- apply(contr, 2, mean) * 100
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
    function(object)
{
    out <- lapply(object, function(x) t(x[x$cum <= 70 ,"cum", drop = FALSE]))
    print(out)
    invisible(object)
}

`summary.simper` <-
    function(object)
{
    class(object) <- "summary.simper"
    object
}

#' Similarity Percentages
#'
#' Discriminating species between two groups using Bray-Curtis dissimilarities.
#'
#' @param comm Community data
#' @param group Vector assigning the groups
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
#' @references Clarke, K.R. 1993. Non-parametric multivariate analyses of changes in community structure. \emph{Austral Ecology}, 18, 117–143.
#' @export
#'
#' @keywords multivariate
#'
#' @examples
#' require(vegan)
#' data(dune)
#' data(dune.env)
#' with(dune.env, simper(dune, Management))

simper <- function(comm, group)
{
    comp <- t(combn(unique(as.character(group)), 2))
    outlist <- NULL
    for (i in 1:nrow(comp)) {
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
                    a <- rowSums(me)
                    c <- md / a
                    contr[(j-1)*n.a+k, ] <- md / a	
                }
            }
        }
        av.contr <- apply(contr, 2, mean) * 100
        ov.av.dis <- sum(av.contr)
        sdi <- apply(contr, 2, sd)
        sdi.av <- sdi / av.contr
        av.a <- colMeans(group.a)
        av.b <- colMeans(group.b) 
        dat <- data.frame(av.contr, sdi, sdi.av, av.a, av.b)
        dat <- dat[order(dat$av.contr, decreasing = TRUE), ]
        cum <-  cumsum(dat$av.contr / ov.av.dis) * 100
        out <-  data.frame(dat, cum)
        names(out) <- c("contr", "sd", "contr/sd", paste("av_", comp[i, 1], sep = ""), paste("av_", comp[i, 2], sep = ""), "cum")
        outlist[[paste(comp[i,1], "_", comp[i,2], sep = "")]] <- out
    }
    outlist
}

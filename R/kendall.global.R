`kendall.global` <-
	function(Y, group, nperm=999, mult="holm")
{
### Function to test the overall significance of the Kendall coefficient of  
### concordance for a single group or several groups of judges (e.g. species)
###
### copyleft - Guillaume Blanchet and Pierre Legendre, October 2008
################################################################################
	
    mult <- match.arg(mult, c("sidak", "holm", "bonferroni"))
	
    a <- system.time({

        ##CC# Make sure Y is a matrix and find number of rows and columns of Y
	Y <- as.matrix(Y)
	n <- nrow(Y)
	p <- ncol(Y)
	if(p < 2) stop("There is a single variable in the data matrix")

        ##CC# Transform the species abundances to ranks, by column
	R <- apply(Y,2,rank)

	if(missing(group)) group <- rep(1,p)
	if(length(group) != p){
            stop("The number of species in the vector differs from the total number of species")
	}
	group <- as.factor(group)
	gr.lev <- levels(group)
	ngr <- nlevels(group)
		
	gr <- as.list(1:ngr)
		
	n.per.gr <- vector(length=ngr)
	for(i in 1:ngr){
            gr[[i]] <- which(group==gr.lev[i])
            n.per.gr[i] <- length(gr[[i]])
            ## Vector containing the number of species per group
	}

        ## Initialise the vectors of results
	W.gr <- vector("numeric",ngr)
	Chi2.gr <- vector("numeric",ngr)
	P.Chi2.gr <- vector("numeric",ngr)
	vec <- NA

	for(i in 1:ngr){
            p.i <- n.per.gr[i]
            if(p.i < 2) stop("There is a single variable in group ",gr.lev[i])
            ##CC# Correction factors for tied ranks (eq. 3.3)
            t.ranks <- apply(R[,gr[[i]]], 2, function(x) summary(as.factor(x)))
            TT <- sum(unlist(lapply(t.ranks, function(x) sum((x^3)-x))))
	
            ##CC# Calculate the Sum of squares of the uncentred ranks (S) (eq. 3.1)
            S <- sum(apply(R[,gr[[i]]], 1, sum)^2)
	
            ##CC# Calculate Kendall's W (eq. 3.2)
            W.gr[i] <- ((12*S)-(3*(p.i^2)*n*(n+1)^2))/(((p.i^2)*((n^3)-n))-(p.i*TT))
	
            ##CC# Calculate Friedman's Chi-square (eq. 3.4)
            Chi2.gr[i] <- p.i*(n-1)*W.gr[i]

            counter <- 1
            for(j in 1:nperm) { # Each species is permuted independently
                R.perm <- apply(R[,gr[[i]]], 2, sample)
                S.perm <- sum(apply(R.perm, 1, sum)^2)
                W.perm <- ((12*S.perm)-(3*(p^2)*n*(n+1)^2))/(((p^2)*((n^3)-n))-(p*TT))
                Chi2.perm <- p*(n-1)*W.perm
                if(Chi2.perm >= Chi2.gr[i]) counter <- counter+1
            }
            P.Chi2.gr[i] <- counter/(nperm+1)
            vec <- c(vec, P.Chi2.gr[i])
	}
        ## Correction to P-values for multiple testing
        ## Write all P-values to vector 'vec'
        if(ngr > 1) {
            vec <- vec[-1]
            if(length(vec) != ngr) stop("Error in putting together vector 'vec'")

            if(mult == "sidak") {
                vec.corr = NA
                for(i in 1:ngr) vec.corr = c(vec.corr, (1-(1-vec[i])^ngr))
                vec.corr <- vec.corr[-1]
            }
            if(mult == "holm") vec.corr <- p.adjust(vec, method="holm")
            if(mult == "bonferroni") vec.corr <- p.adjust(vec, method="bonferroni")

            vec.gr <- vector("list",ngr)
            end <- 0
            for(i in 1:ngr){
                beg <- end+1
                end <- end + n.per.gr[i]
                vec.gr[[i]] <- vec.corr[beg:end]
            }
        }
        ## Create a data frame containing the results
	if(ngr == 1) {
            table <- rbind(W.gr, Chi2.gr, P.Chi2.gr)
            colnames(table) <- colnames(table,do.NULL = FALSE, prefix = "Group.")
            rownames(table) <- c("W", "Chi2", "Prob")
        }  else {
            table <- rbind(W.gr, Chi2.gr, P.Chi2.gr, vec.corr)
            colnames(table) <- colnames(table,do.NULL = FALSE, prefix = "Group.")
            rownames(table) <- c("W", "Chi2", "Prob", "Corrected prob")		
        }
                                        #
    })
    a[3] <- sprintf("%2f",a[3])
    cat("Time to compute global test(s) =",a[3]," sec",'\n')
                                        #
    if(ngr == 1) {
        out <- list(Concordance_analysis=table)
    }  else {
        out <- list(Concordance_analysis=table, Correction.type=mult)
    }
                                        #
    class(out) <- "kendall.global"
    out
}

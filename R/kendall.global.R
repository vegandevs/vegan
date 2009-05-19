`kendall.global` <-
	function(Y, group, nperm=999, mult="holm")
{
### Function to test the overall significance of the Kendall coefficient of  
### concordance for a single group or several groups of judges (e.g. species)
###
### copyleft - Guillaume Blanchet and Pierre Legendre, October 2008
################################################################################
	
    mult <- match.arg(mult, c("sidak", p.adjust.methods))
	
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
    F.gr <- vector("numeric",ngr)
    prob.F.gr <- vector("numeric",ngr)
    Chi2.gr <- vector("numeric",ngr)
    prob.perm.gr <- vector("numeric",ngr)
    vec <- NA

    for(i in 1:ngr){
        p.i <- n.per.gr[i]
        if(p.i < 2) stop("There is a single variable in group ",gr.lev[i])
        ##CC# Correction factors for tied ranks (eq. 3.3)
        t.ranks <- apply(R[,gr[[i]]], 2,
                         function(x) summary(as.factor(x), maxsum=n))
        T. <- sum(unlist(lapply(t.ranks, function(x) sum((x^3)-x))))
	
        ##CC# Compute the Sum of squares of the uncentred ranks (S) (eq. 3.1)
        S <- sum(apply(R[,gr[[i]]], 1, sum)^2)
	
        ##CC# Compute Kendall's W (eq. 3.2)
        W.gr[i] <- ((12*S)-(3*(p.i^2)*n*(n+1)^2))/(((p.i^2)*((n^3)-n))-(p.i*T.))
		
        ##C# Compute Fisher F statistic and associated probability
        F.gr[i] <- (p.i-1)*W.gr[i]/(1-W.gr[i])
        nu1 <- n-1-(2/p.i)
        nu2 <- nu1*(p.i-1)
        prob.F.gr[i] <- pf(F.gr[i], nu1, nu2, lower.tail=FALSE)
	
        ##CC# Calculate Friedman's Chi-square (eq. 3.4)
        Chi2.gr[i] <- p.i*(n-1)*W.gr[i]

        counter <- 1
        for(j in 1:nperm) {   # Each species is permuted independently
            R.perm <- apply(R[,gr[[i]]], 2, sample)
            S.perm <- sum(apply(R.perm, 1, sum)^2)
            W.perm <- ((12*S.perm)-(3*(p.i^2)*n*(n+1)^2))/(((p.i^2)*((n^3)-n))-(p.i*T.))
            Chi2.perm <- p.i*(n-1)*W.perm
            if(Chi2.perm >= Chi2.gr[i]) counter <- counter+1
        }
        prob.perm.gr[i] <- counter/(nperm+1)
    }
    ## Correction to P-values for multiple testing
    if(ngr > 1) {
        if(mult == "sidak") {
            perm.corr <- NA
            for(i in 1:ngr) perm.corr = c(perm.corr, (1-(1-prob.perm.gr[i])^ngr))
            perm.corr <- perm.corr[-1]
                                        #
            prob.F.corr <- NA
            for(i in 1:ngr) prob.F.corr = c(prob.F.corr, (1-(1-prob.F.gr[i])^ngr))
            prob.F.corr <- prob.F.corr[-1]
        } else {
            perm.corr <- p.adjust(prob.perm.gr, method=mult)
            prob.F.corr <- p.adjust(prob.F.gr, method=mult)
        }
    }
    ## Create a data frame containing the results
    if(ngr == 1) {
        table <- rbind(W.gr, F.gr, prob.F.gr, Chi2.gr, prob.perm.gr)
        colnames(table) <- colnames(table,do.NULL = FALSE, prefix = "Group.")
        rownames(table) <- c("W", "F", "Prob.F", "Chi2", "Prob.perm")
    }  else {
        table <- rbind(W.gr, F.gr, prob.F.gr, prob.F.corr, Chi2.gr, prob.perm.gr, perm.corr)
        colnames(table) <- colnames(table,do.NULL = FALSE, prefix = "Group.")
        rownames(table) <- c("W", "F", "Prob.F", "Corrected prob.F", "Chi2", "Prob.perm", "Corrected prob.perm")		
    }
    if(ngr == 1) {
        out <- list(Concordance_analysis=table)
    }  else {
        out <- list(Concordance_analysis=table, Correction.type=mult)
    }
                                        #
    class(out) <- "kendall.global"
    out
}

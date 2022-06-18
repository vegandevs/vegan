`kendall.post` <-
	function(Y, group, nperm=999, mult="holm")
{
### Function to carry out a posteriori tests on individual judges (e.g. species)
### for a single group or several groups of judges.
###
### copyleft - Guillaume Blanchet and Pierre Legendre, October 2008
################################################################################

    mult <- match.arg(mult, c("sidak", p.adjust.methods))

    ##CC# Make sure Y is a matrix and find number of rows and columns of Y
    Y <- as.matrix(Y)
    p <- ncol(Y)
    if(p < 2) stop("there is only one variable in the data matrix")

    ##CC# Transform the species abundances to ranks, by column
    R <- apply(Y,2,rank)

    if(missing(group)) group <- rep(1,p)
    if(length(group) != p){
        stop("the number of species in the vector differs from the total number of species")
    }

    ##CC# Separate tests for the variables in each group
    group <- as.factor(group)
    gr.lev <- levels(group)
    ngr <- nlevels(group)

    gr <- as.list(1:ngr)

    n.per.gr <- vector(length=ngr)
    for(i in 1:ngr){
        gr[[i]] <- which(group==gr.lev[i])
        n.per.gr[i] <- length(gr[[i]])  # Vector with the number of
                                        # variables per group
    }

    ##===============================
    ##CC# start permutation procedure
    ##===============================
    ##CC# Initialize the counters
    counter <- as.list(1:ngr)
    for(i in 1:ngr){
        counter[[i]] <- vector(length = n.per.gr[i])
    }
    W.gr <- vector("list",ngr)
    if(ngr > 1) spear.gr <- vector("list",ngr)

    for(i in 1:ngr){
        p.i <- n.per.gr[i]
        if(p.i < 2) stop(gettextf("there is only one variable in group %d",
                                  gr.lev[i]))
                                        #CC# Extract variables part of
                                        #group i
        R.gr <- R[,gr[[i]]]     # Table with species of group 'i' only
                                        #as columns CC# Calculate the
                                        #mean of the Spearman
                                        #correlations for each species
                                        #in a group
        spear.mat <- cor(R.gr)
        diag(spear.mat) <- NA
        spear.mean <- apply(spear.mat, 1, mean, na.rm=TRUE)
                                        #CC# Calculate Kendall's W for each variable
        W.var <- ((p.i-1)*spear.mean+1)/p.i

                                        #for(j in 1:n.per.gr[i]){
        for(j in 1:p.i){                # Test each species in turn
                                        #CC# Create a new R where the
                                        #permuted variable has been
                                        #removed
            R.gr.mod <- R.gr[,-j]
            ##CC# Set counter
            counter[[i]][j] <- 1
            ##CC# The actual permutation procedure
            for(k in 1:nperm){
                var.perm <- sample(R.gr[,j])
                spear.mat.perm <- cor(cbind(R.gr.mod, var.perm))
                diag(spear.mat.perm) <- NA
                spear.mean.j.perm <- mean(spear.mat.perm[p.i,], na.rm=TRUE)
                W.perm <- ((p.i-1)*spear.mean.j.perm+1)/p.i
                if(W.perm >= W.var[j]) counter[[i]][j] <- counter[[i]][j]+1
            }
        }
        W.gr[[i]] <- W.var
        if(ngr > 1) spear.gr[[i]] <- spear.mean
    }
    ##CC# Calculate P-values
    for(i in 1:ngr) {
        counter[[i]] <- counter[[i]]/(nperm+1)
    }

    ## Correction to P-values for multiple testing
    ## Write all P-values to a long vector 'vec'
    vec <- counter[[1]]
    if(ngr > 1) {
        for(i in 2:ngr) vec = c(vec, counter[[i]])
    }
    if(length(vec) != p) stop("error in putting together vector 'vec'")

    if(mult == "sidak") {
        vec.corr = NA
        for(i in 1:p) vec.corr = c(vec.corr, (1-(1-vec[i])^p))
        vec.corr <- vec.corr[-1]
    } else {
        vec.corr <- p.adjust(vec, method=mult)
    }

    if(ngr > 1) {
        vec.gr <- vector("list",ngr)
        end <- 0
        for(i in 1:ngr){
            beg <- end+1
            end <- end + n.per.gr[i]
            vec.gr[[i]] <- vec.corr[beg:end]
        }
    }
    ## Create data frames containing the results
    if(ngr == 1) {
        table <- rbind(spear.mean, W.var, counter[[1]], vec.corr)
        rownames(table) <- c("Spearman.mean", "W.per.species", "Prob", "Corrected prob")
    } else {
        table <- as.list(1:ngr)
        for(i in 1:ngr) {
            table[[i]] <- rbind(spear.gr[[i]], W.gr[[i]], counter[[i]], vec.gr[[i]])
            rownames(table[[i]]) <- c("Spearman.mean", "W.per.species", "Prob", "Corrected prob")
            ## PL: Next line had been lost
            colnames(table[[i]]) <- colnames(table[[i]], do.NULL = FALSE,
                                             prefix = "Spec")
        }
    }
    if(ngr == 1) {
        out <- list(A_posteriori_tests=table, Correction.type=mult)
    } else {
        out <- list(A_posteriori_tests_Group=table, Correction.type=mult)
    }
                                        #
    class(out) <- "kendall.post"
    out
}

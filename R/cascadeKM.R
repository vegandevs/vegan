"cascadeKM" <-
function(data, inf.gr, sup.gr, iter = 100, criterion="calinski")
{
### DESCRIPTION

### This function use the 'kmeans' function of the 'stats' package to create 
### a cascade of partitions from K = nb_inf_gr to K = nb_sup_gr

### INPUT				
###
### data			The data matrix; the objects are the rows
### nb_inf_gr		Number of groups (K) for the first partition (min)
### nb_sup_gr	 	Number of groups (K) for the last partition (max)
### iter			The number of random starting configurations for each value of K
### criterion		The criterion that will be used to select the best
###				partition. See the 'clustIndex' function in PACKAGE = cclust

### OUTPUT
###
### The same as in the kmeans packages

### EXAMPLE
###
### 	result <- cascadeKM(donnee, 2, 30, iter = 50, criterion = 'calinski') 
###
### 	data = data table
### 	2 = lowest number of groups for K-means
### 	30 = highest number of groups for K-means
### 	iter = 50: start kmeans 50 times using different random configurations
### 	criterion = 'calinski': the Calinski-Harabasz (1974) criterion to determine 
###      the best value of K for the data set. 'Best' is in the least-squares sense.
###

### Main function
    data <- as.matrix(data)
    SCE<-list()
    resultat<-list()
    index<-list()
    if(!is.null(nrow(data))){
        partition <- matrix(NA, nrow(data), sup.gr - inf.gr + 1)
    }else{
        partition <- matrix(NA, length(data), sup.gr - inf.gr + 1)
    }
    results <- matrix(NA, 2, sup.gr - inf.gr + 1)
    size <- matrix(NA, sup.gr, sup.gr - inf.gr + 1)
    ## Pour tous les nombres de groupes voulus
    h <- 1
    for(ii in inf.gr:sup.gr)
    {
        ## Initialization
        ## set.seed(ii)  
        ## Set.seed ˆ ŽtŽ enlevŽ car il rend instable la fonction kmeans
        j <- ii - inf.gr + 1
        tmp <- kmeans(data, ii, iter.max = 50, nstart=iter)
        size[1:ii,h] <- tmp$size
        h <- h+1
        partition[,j] <- tmp$cluster
        ## Compute SSE statistic
        results[1,j] <- sum(tmp$withinss)
        ## Compute stopping criterion
        results[2,j] <- cIndexKM(tmp,data, index = tolower(criterion))
    }
    colnames(partition) <- paste(inf.gr:sup.gr, "groups")
    tmp <- rownames(data)
    if(is.null(tmp)){
        r.name <- c(1:nrow(partition))
    }else{
        r.name <- tmp
    }
    rownames(partition) <- r.name
	
    colnames(results) <- paste(inf.gr:sup.gr, "groups")
    rownames(results)<-c("SSE", criterion)

    colnames(size) <- paste(inf.gr:sup.gr, "groups")
    rownames(size) <- paste("Group", 1:sup.gr)
	
    tout<-list(partition=partition, results=results, criterion=criterion, size=size)
    class(tout) <- "cascadeKM"
    tout
}

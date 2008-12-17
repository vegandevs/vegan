### NODF metric of nestedness.
### Code submitted as an R-Forge feature request #265 by Gustavo
### Carvalho.
nestednodf <- function(comm)
{
    ## Transform the community matrix to presence/ausence  
    comm <- ifelse(comm > 0, 1, 0)
    ## Order rows and columns
    comm <- comm[order(rowSums(comm), decreasing=TRUE),
                 order(colSums(comm), decreasing=TRUE)]  
    ## Two matrices with all possible rows and columns combinations
    ## used to calculate the paired overlap
    row.combinations <- t(combn(1:dim(comm)[1],2))
    col.combinations <- t(combn(1:dim(comm)[2],2))
    ## Saving a bit of cpu time
    row.sums <- rowSums(comm)
    col.sums <- colSums(comm)
    dimensions <- dim(comm)  
    fill <- sum(comm)/length(comm)

    ## N.paired for columns
    combcol <- function(x) {
        N.paired <- 0
        ## if MTi > MTj
        if (diff(col.sums[x]) < 0 && !any(col.sums[x] == 0)) {     
            paired.overlap <- sum(rowSums(comm[,x]) == 2) / col.sums[x[2]]
            N.paired <- paired.overlap
        }
        return(N.paired)
    }

    ## N.paired for rows
    combrow <- function(x) {
        N.paired <- 0
        ## If MTk > MTl
        if (diff(row.sums[x]) < 0 && !any(row.sums[x] == 0)) {
            paired.overlap <- sum(colSums(comm[x,]) == 2) / row.sums[x[2]]
            N.paired <- paired.overlap
        }
        return(N.paired)
    }

    ## N.paired for all columns and rows
    N.paired.columns <- apply(col.combinations, 1, combcol)
    N.paired.rows <- apply(row.combinations, 1, combrow) 

    ## NODF, N.rows, and N.columns
    NODF <- (sum(c(N.paired.rows, N.paired.columns)) * 100)/
        ((dimensions[2] * (dimensions[2] - 1) / 2) +
         (dimensions[1] * (dimensions[1] - 1) / 2))
    N.rows <- mean(N.paired.rows) * 100  
    N.columns <- mean(N.paired.columns) * 100
  
    ## Returned list (changed to make it more similar to what
    ## nestedtemp returns).
    out <- list(comm = comm, fill = fill, N.rows = N.rows,
                N.columns = N.columns, statistic = NODF)
    return(out)
}

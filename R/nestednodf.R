`nestednodf` <- 
    function(comm, order = TRUE)
{
    comm <- ifelse(comm > 0, 1, 0)
    ## Order rows and columns
    if (order)
        comm <- comm[order(rowSums(comm), decreasing=TRUE),
                     order(colSums(comm), decreasing=TRUE)]    
    dimensions <- dim(comm)
    fill <- sum(comm)/length(comm)
    N.paired <- 0
    ## Function to be applied to each combination of rows and columns
    comb <- function(x, rows) {
        if (identical(rows,TRUE)) {
            comb.first <- comm[x[1],]
            comb.second <- comm[x[2],]
        }
        else {
            comb.first <- comm[,x[1]]
            comb.second <- comm[,x[2]]
        }
        ## if MTi > MTj
        if (sum(comb.first) > sum(comb.second) && sum(comb.second) > 0) {
            paired.overlap <- sum((comb.first + comb.second) == 2) /
                sum(comb.second)
            N.paired <- paired.overlap
        }
        return(N.paired)
    }
    ## N.paired for all combinations of columns and rows
    N.paired.rows <- combn(1:dimensions[1],2, comb, rows=TRUE)
    N.paired.columns <- combn(1:dimensions[2],2, comb, rows=FALSE)
    ## Index calculations
    N.columns <- mean(N.paired.columns) * 100
    N.rows <- mean(N.paired.rows) * 100  
    NODF <- (sum(c(N.paired.rows, N.paired.columns)) * 100) /
        ((dimensions[2] * (dimensions[2] - 1) / 2) + 
         (dimensions[1] * (dimensions[1] - 1) / 2))
    ## Returned list
    out <- list(comm = comm, fill = fill, 
                statistic=c("N.columns" = N.columns, "N.rows" = N.rows, "NODF" = NODF))
    class(out) <- "nestednodf"
    return(out)
}

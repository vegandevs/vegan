`nestednodf` <- 
    function(comm, order = TRUE, weighted = FALSE)
{
    ## Keep it consistent with similar functions
    if (!weighted) {
        comm <- ifelse(comm > 0, 1, 0)
    }
    if (order) {
        comm <- comm[order(rowSums(comm), decreasing=TRUE),
                     order(colSums(comm), decreasing=TRUE)]    
    }
    fill <- comm > 0
    rfill <- rowSums(fill)
    cfill <- colSums(fill)
    nr <- NROW(comm)
    nc <- NCOL(comm)
    fill <- sum(rfill)/length(comm)
    ## Initializing the vectors that will hold nodf values
    ## for each combination of rows and columns
    N.paired.rows <- numeric(nr * (nr - 1) / 2)
    N.paired.cols <- numeric(nc * (nc - 1) / 2)
    counter <- 0
    ## Nested loops to get every combination of rows/columns
    ## Row i
    for (i in 1:(nr-1)) {
        first <- comm[i, ]
        ## Row j
        for(j in (i + 1):nr) {
            ## Counting the number of nodfs calculated so far
            counter <- counter + 1
            if (rfill[i] <= rfill[j] || any(rfill[c(i, j)] == 0))
                next
            if (weighted) {
                second <- comm[j, ]
                N.paired.rows[counter] <- sum(first - second > 0 & second > 0) /
                    sum(second > 0)
            }
            else {
                N.paired.rows[counter] <- sum(first + comm[j, ] == 2) /
                    rfill[j]
            }
        }
    }
    ## Reseting the counter
    counter <- 0
    ## Column i
    for (i in 1:(nc-1)) {
        first <- comm[, i]
        ## Column j
        for(j in (i + 1):nc) {
            counter <- counter + 1
            if (cfill[i] <= cfill[j] || any(cfill[c(i, j)] == 0))
                next
            if (weighted) {
                second <- comm[, j]
                N.paired.cols[counter] <- sum(first - second > 0 & second > 0) /
                    sum(second > 0)
            }
            else {
                N.paired.cols[counter] <- sum(first + comm[, j] == 2) /
                    cfill[j]
            }
        }
    }
    ## Calculating the means and the NODF for the whole matrix
    N.columns <- mean(N.paired.cols) * 100
    N.rows <- mean(N.paired.rows) * 100
    NODF <- (sum(c(N.paired.rows, N.paired.cols)) * 100)/
        ((nc * (nc - 1)/2) + (nr * (nr - 1)/2))
    out <- list(comm = comm, fill = fill, 
                statistic=c("N.columns" = N.columns, "N.rows" = N.rows,
                "NODF" = NODF))
    class(out) <- "nestednodf"
    return(out)
}

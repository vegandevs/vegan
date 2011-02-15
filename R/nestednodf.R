`nestednodf` <- 
    function(comm, order = TRUE, weighted = FALSE) 
{
    bin.comm <- ifelse(comm > 0, 1, 0)
    rfill <- rowSums(bin.comm)
    cfill <- colSums(bin.comm)
    if (!weighted)
        comm <- bin.comm
    if (order) {
        if (weighted) {
            rgrad <- rowSums(comm)
            cgrad <- colSums(comm)
            rorder <- order(rfill, rgrad, decreasing = TRUE)
            corder <- order(cfill, cgrad, decreasing = TRUE)
        } else {
            rorder <- order(rfill, decreasing = TRUE)
            corder <- order(cfill, decreasing = TRUE)
        }
        comm <- comm[rorder, corder]
        rfill <- rfill[rorder]
        cfill <- cfill[corder]
    }
    nr <- NROW(comm)
    nc <- NCOL(comm)
    fill <- sum(rfill)/length(comm)
    N.paired.rows <- numeric(nr * (nr - 1)/2)
    N.paired.cols <- numeric(nc * (nc - 1)/2)
    counter <- 0
    for (i in 1:(nr - 1)) {
        first <- comm[i, ]
        for (j in (i + 1):nr) {
            counter <- counter + 1
            if (rfill[i] <= rfill[j] || any(rfill[c(i, j)] == 0)) 
                next
            if (weighted) {
                second <- comm[j, ]
                N.paired.rows[counter] <-
                    sum(first - second > 0 & second > 0)/sum(second > 0)
            }
            else {
                N.paired.rows[counter] <-
                    sum(first + comm[j, ] == 2)/rfill[j]
            }
        }
    }
    counter <- 0
    for (i in 1:(nc - 1)) {
        first <- comm[, i]
        for (j in (i + 1):nc) {
            counter <- counter + 1
            if (cfill[i] <= cfill[j] || any(cfill[c(i, j)] == 0)) 
                next
            if (weighted) {
                second <- comm[, j]
                N.paired.cols[counter] <-
                    sum(first - second > 0 & second > 0)/sum(second > 0)
            }
            else {
                N.paired.cols[counter] <-
                    sum(first + comm[, j] == 2)/cfill[j]
            }
        }
    }
    N.columns <- mean(N.paired.cols) * 100
    N.rows <- mean(N.paired.rows) * 100
    NODF <- (sum(c(N.paired.rows, N.paired.cols)) * 100)/
        ((nc * (nc - 1)/2) + (nr * (nr - 1)/2))
    out <- list(comm = comm, fill = fill,
                statistic = c(N.columns = N.columns, N.rows = N.rows, NODF = NODF))
    class(out) <- "nestednodf"
    return(out)
}

`read.cep` <-
    function (file,  positive = TRUE, trace = FALSE)
{
    ## skip first line and get the format
    cep <- readLines(file)
    i <- 2
    fmt <- substr(cep[i], 1, 60)
    fmt <- toupper(fmt)
    fmt <- gsub(" ", "", fmt)
    ## get the number data entries per record
    nrecord <- as.numeric(substr(cep[i], 61,80))
    if (is.na(nrecord)) {
        i <- i+1
        nrecord <- as.numeric(cep[i])
    }
    ## process format
    fmt <- gsub(paste0(nrecord, "\\("), "", fmt)
    fmt <- gsub("\\(","",fmt)
    fmt <- gsub("\\)","",fmt)
    fmt <- strsplit(fmt, ",")[[1]]
    len <- length(fmt)
    fmt <- c(fmt[seq_len(len-2)], rep(fmt[c(len-1,len)], nrecord))
    ## vectors to store results (with safe size)
    nlines <- length(cep)-i
    siteid <- numeric(nlines * nrecord)
    specid <- numeric(nlines * nrecord)
    abund <- numeric(nlines * nrecord)
    ids <- seq(2, by=2, len=nrecord)
    id <- 0
    ## read until there an empty siteid
    repeat {
        i <- i+1
        x <- drop(as.matrix(read.fortran(textConnection(cep[i]), fmt)))
        if(x[1] <= 0) break
        for(j in ids) {
            if(!is.na(x[j])) {
                id <- id+1
                siteid[id] <- x[1]
                specid[id] <- x[j]
                abund[id] <- x[j+1]
            }
            else
                break
        }
    }
    ## max identifiers
    nsp <- max(specid)
    nst <- max(siteid)
    ## make as a matrix
    out <- matrix(0, nst, nsp)
    for(j in seq_len(id))
        out[siteid[j], specid[j]] <- abund[j]
    out
}

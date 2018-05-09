### Reads condensed CEP and similar CANOCO formatted data

### This function was originally based on Fortran code to interpret
### the format and read the data, but Fortran I/O is no longer allowed
### in CRAN packages. The original Fortran version was made available
### in github package vegandevs/cepreaded. This function uses
### utils::read.fortran. For that function we must translate the
### Fortran format in form understood by read.fortran(). This may
### fail, and read.fortran() can also fail, in particular wih
### interpreting the length of decimal part in F format.

### The original Fortran (and cepreader) can also interpret open
### (ordinary) Fortran format and Canoco FREE format, but this
### function currently only reads condensed data.
`read.cep` <-
    function (file,  positive = TRUE)
{
    ## Read all file first, interpret contents later
    cep <- readLines(file)
    ## skip first line and get the format
    i <- 2
    fmt <- substr(cep[i], 1, 60) # CEP sets fmt in cols 1-60
    fmt <- toupper(fmt)
    fmt <- gsub(" ", "", fmt)
    ## get the number of data entries per record
    nrecord <- as.numeric(substr(cep[i], 61,80))
    if (is.na(nrecord)) {
        i <- i+1
        nrecord <- as.numeric(cep[i])
    }
    ## check format: there should be to I-elements (site id, species
    ## id), and there should be two opening "("
    fmt1 <- strsplit(fmt, NULL)[[1]]
    if (sum(fmt1 == "I") != 2 || (nrecord > 1 && sum(fmt1 == "(") != 2))
        stop(gettextf("format %s does not look correct for condensed data",
                      fmt))
    ## process format: basically the format should have elements for
    ## (INT, n(INT, REAL)). read.fortran() does not understand
    ## multiplier 'n', but we need to rep((INT,REAL), n) in the
    ## fmt vector.
    fmt <- gsub(paste0(nrecord, "\\("), ";", fmt) # separate with ;
    fmt <- gsub("\\(","",fmt)
    fmt <- gsub("\\)","",fmt)
    ## number of decimals: there should be one and only one Fa.b
    ## format, and we need 'b'
    ndec <- as.numeric(strsplit(fmt, "\\.")[[1]][2])
    ## now split format for plotid and nrecord couplets
    fmt <- strsplit(fmt, ";")[[1]]
    fmt <- c(strsplit(fmt[1], ",")[[1]], rep(strsplit(fmt[2], ",")[[1]],
                                             nrecord))
    if (any(is.na(fmt)))
        fmt <- fmt[!is.na(fmt)]
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
        if(is.na(x[1]) || x[1] <= 0) break
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
    ## check there are no duplicate entries: only last one would be
    ## used, and this causes an error (and this has happened in
    ## literature)
    if (any(dups <- duplicated(cbind(siteid, specid))[seq_len(id)]))
        stop("you have duplicated data entries: ",
             paste(siteid[seq_len(id)][dups], specid[seq_len(id)][dups],
                   collapse = ", "))
    ## max identifiers
    nsp <- max(specid)
    nst <- max(siteid)
    ## read dimnames
    i <- i+1
    nomina <- read.fwf(textConnection(cep[i:length(cep)]), rep(8, 10),
                       as.is=TRUE)
    nomina <- gsub(" ", "", as.vector(t(nomina)))
    spnam <- make.names(nomina[seq_len(nsp)], unique = TRUE)
    nst0 <- ceiling(nsp/10) * 10
    stnam <- make.names(nomina[seq_len(nst) + nst0], unique = TRUE)
    ## utils::read.fortran divides with 10^ndec of F format even when
    ## there is an explicit decimal point: undo if this seems to have
    ## happened
    if (ndec > 0 && min(abund[1:id]) <= 10^(-ndec))
        abund <- abund * 10^ndec
    ## make as a matrix
    out <- matrix(0, nst, nsp)
    for(j in seq_len(id))
        out[siteid[j], specid[j]] <- abund[j]
    dimnames(out) <- list(stnam, spnam)
    if (positive) {
        rs <- rowSums(out)
        cs <- colSums(out)
        if (any(cs <= 0) || any(rs <= 0))
            out <- out[rs > 0, cs > 0]
    }
    as.data.frame(out)
}

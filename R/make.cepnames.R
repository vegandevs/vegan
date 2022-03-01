`make.cepnames` <-
    function (names, seconditem = FALSE)
{
    ## make valid names
    names <- make.names(names, unique = FALSE, allow_ = FALSE)
    ## remove trailing and duplicated dots
    names <- gsub("\\.[\\.]+", ".", names)
    names <- gsub("\\.$", "", names)
    ## split by dots and take 4 letters of each element (if several)
    names <- lapply(strsplit(names, "\\."), function(x) if (length(x) > 1)
                    substring(x, 1, 4) else x )
    ## Take first and last element or 8 characters if only one element
    names <- unlist(lapply(names, function(x) if (length(x) > 1)
                           paste(x[c(1, if(seconditem) 2 else length(x))], collapse = "")
                           else x))
    names <- abbreviate(names, 8)
    ## Final clean-up
    make.names(names, unique = TRUE)
}

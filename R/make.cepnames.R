"make.cepnames" <-
    function (names) 
{
    names <- make.names(names, unique = FALSE)
    names <- lapply(strsplit(names, "\\."), function(x) if (length(x) > 1) 
                    substring(x, 1, 4) else x )
    names <- unlist(lapply(names, function(x) if (length(x) > 1)
                           paste(x[c(1, length(x))], collapse = "")
                           else x))
    names <- abbreviate(names, 8)
    names <- make.names(names, unique = TRUE)
    names
}

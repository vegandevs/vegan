`make.cepnames` <-
    function (names, seconditem = FALSE, uniqgenera = FALSE)
{
    ## do not split by hyphens, but collapse hyphened names
    names <- gsub("-", "", names)
    ## make valid names
    names <- make.names(names, unique = FALSE, allow_ = FALSE)
    ## remove trailing and duplicated dots
    names <- gsub("\\.[\\.]+", ".", names)
    names <- gsub("\\.$", "", names)
    ## split by dots and get genus and epithet
    names <- strsplit(names, ".", fixed = TRUE)
    gen <- sapply(names, function(x) x[1])
    epi <- sapply(names,
                  function(x) {if (seconditem) x[2]
                               else if (length(x) > 1) x[length(x)] else ""})
    ## strict=TRUE always takes 4 even if this generates duplicates
    gen <- ifelse(epi != "",
                  abbreviate(gen, 4, use.classes = FALSE, strict = !uniqgenera),
                  gen)
    if (uniqgenera)
        gen <- abbreviate(gen, 4, use.classes = TRUE)
    names <- abbreviate(paste0(gen, epi), 8, use.classes = FALSE)
    ## try to remove wovels if names > 8
    names <- abbreviate(names, 8, use.classes = TRUE, named = FALSE)
    names
}

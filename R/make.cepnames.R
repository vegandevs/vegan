`make.cepnames` <-
    function (names, minlengths = c(4,4), seconditem = FALSE,
              uniqgenera = FALSE, named = FALSE, method)
{
    orignames <- names
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
    ## check duplicates of input names
    if (any(dup <- duplicated(orignames)))
        warning("some input names were duplicated: ",
                paste(orignames[dup], collapse = ", "))
    else if (any(dup <- duplicated(paste(gen, epi))))
        warning("some inferred binomial names were duplicated: ",
                paste(paste(gen, epi)[dup], collapse = ", "))
    ## strict=TRUE always takes given minlength even if these are duplicates
    glen <- minlengths[1]
    nmlen <- sum(minlengths)
    if (missing(method))
        method <- "left.kept"
    gen <- ifelse(epi != "",
                  abbreviate(abbreviate(gen, glen, use.classes = FALSE,
                                        strict = !uniqgenera),
                             glen, use.classes = TRUE, method = method),
                  gen)
    names <- abbreviate(paste0(gen, epi), nmlen, use.classes = FALSE)
    ## try to remove wovels if names > nmlen
    names <- abbreviate(names, nmlen, use.classes = TRUE, method = method,
                        named = FALSE)
    if (named)
        names(names) <- orignames
    ## break duplicates by adding number to the name
    if (any(dup <- duplicated(names))) {
        names <- make.unique(names, sep = "")
        warning("these cepnames were made unique by appending a number:\n",
                paste(paste0(names[dup], " (", orignames[dup], ")"),
                      collapse = ", "))
    }
    names
}

### summary allocates shared components equally to explanatory sets

`summary.varpart` <-
    function(object, ...)
{
    nsets <- object$part$nsets
    nfract <- 2^nsets - 1
    ## find fractions each component is made of
    contr <- matrix(0, nfract, nsets,
                    dimnames = list(paste0("[", letters[seq_len(nfract)], "]"),
                                    paste0("X", seq_len(nsets))))
    setfrac <- rownames(object$part$fract)[seq_len(nsets)]
    for (i in seq_len(nsets))
        contr[,i] <- as.numeric(letters[seq_len(nfract)] %in%
                                strsplit(setfrac[i], "")[[1]])
    ## allocate shared fractions equally to all sets
    contr <- sweep(contr, 1, rowSums(contr), "/")
    uniqpart <- object$part$indfract[seq_len(nfract), 3]
    contr <- uniqpart * contr
    ## returned components
    contrpart <- colSums(contr)
    uniqpart <- uniqpart[seq_len(nsets)]
    names(uniqpart) <- names(contrpart)
    out <- list("uniqpart" = uniqpart,
                "contribpart" = contrpart,
                "contributions" = contr,
                "setlabels" = object$tables)
    class(out) <- "summary.varpart"
    out
}

`print.summary.varpart` <-
    function(x, digits = 3, zero.print = "", ...)
{
    ## collect table
    df <- data.frame("Unique" = x$uniqpart, "Contributed" = x$contribpart,
                     "Component" = x$setlabels)
    cat("\nUnique fractions and total with shared fractions equally allocated:\n\n")
    print(df, digits = digits)
    cat("\nContributions of fractions to sets:\n\n")
    print.table(x$contributions, digits = digits, zero.print = zero.print)
    invisible(x)
}

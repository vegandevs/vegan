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
    out <- list("uniqpart" = uniqpart,
                "contribpart" = contrpart,
                "contributions" = contr)
    class(out) <- "summary.varpart"
    out
}

`ordilattice.getEnvfit` <-
    function(formula, object, envfit, choices = 1:3,  ...)
{
    if (!missing(envfit) && !is.null(envfit))
        object <- envfit
    bp <- scores(object, display = "bp", choices = choices)
    cn <- scores(object, display = "cn",  choices = choices)
    bp <- bp[!(rownames(bp) %in% rownames(cn)),, drop=FALSE]
    left <- as.character(formula[[2]])
    right <- formula[[3]]
    if (length(right) == 3)
        right <- right[[2]]
    right <- as.character(right)
    if (all(c(left,right) %in% colnames(bp)))
        bp <- bp[, c(left,right), drop=FALSE]
    else
        bp <- NULL
    if (!is.null(bp) && nrow(bp) == 0)
        bp <- NULL
    if (!is.null(ncol(cn)) && all(c(left,right) %in% colnames(cn)))
        cn <- cn[, c(left,right), drop=FALSE]
    else
        cn <- NULL
    list(arrows = bp, centres = cn)
}

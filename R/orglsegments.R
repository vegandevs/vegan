"orglsegments" <-
    function (object, groups, display = "sites", choices = 1:3,...) 
{
    pts <- scores(object, display = display, choices = choices, ...)
    inds <- names(table(groups))
    for (is in inds) {
        X <- pts[groups == is, , drop = FALSE]
        if (nrow(X) > 1) {
            for (i in 2:nrow(X)) {
                rgl.lines(c(X[i-1,1],X[i,1]), c(X[i-1,2],X[i,2]), 
                          c(X[i-1,3],X[i,3]), ...)
            }
        }
    }
    invisible()
}


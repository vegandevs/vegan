`multipart.formula` <-
    function(formula, data, index=c("renyi", "tsallis"), scales = 1,
             global = FALSE, relative = FALSE, nsimul=99, ...)
{
    ## evaluate formula
    lhs <- formula[[2]]
    if (missing(data))
        data <- parent.frame()
    lhs <- as.matrix(eval(lhs, data))
    formula[[2]] <- NULL
    rhs <- model.frame(formula, data, drop.unused.levels = TRUE)

    ## part check proper design of the model frame
    noint <- attr(attr(attr(rhs, "terms"), "factors"), "dimnames")[[1]]
    int <- attr(attr(attr(rhs, "terms"), "factors"), "dimnames")[[2]]
    if (!identical(noint, int))
        stop("interactions are not allowed in formula")
    if (!all(attr(attr(rhs, "terms"), "dataClasses") == "factor"))
        stop("all right hand side variables in formula must be factors")
    sim <- multipart.default(lhs, rhs, index = index, scales = scales,
                             global = global, relative = relative,
                             nsimul = nsimul, ...)
    call <- match.call()
    call[[1]] <- as.name("multipart")
    attr(sim, "call") <- call
    sim
}

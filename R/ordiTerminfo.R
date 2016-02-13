"ordiTerminfo" <-
    function(d, data)
{
    Terms <- delete.response(d$terms.expand)
    if (length(attr(Terms, "term.labels")) == 0)
        mf <- data.frame(NULL)
    else
        mf <- d$modelframe
    xlev <- .getXlevels(Terms, mf)
    ordered <- sapply(mf, is.ordered)
    list(terms = Terms, xlev = xlev, ordered = ordered)
}

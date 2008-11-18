"ordiTerminfo" <-
    function(d, data)
{
    Terms <- delete.response(d$terms.expand)
    if (length(attr(Terms, "term.labels")) == 0)
        mf <- data.frame(NULL)
    else
        mf <- model.frame(formula(Terms), data)
    xlev <- .getXlevels(Terms, mf)
    ordered <- sapply(mf, is.ordered)
    list(terms = Terms, xlev = xlev, ordered = ordered)
}

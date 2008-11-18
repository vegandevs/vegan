`permControl` <- function(strata = NULL, nperm = 199, complete = FALSE,
                          #type = c("free", "series", "grid", "strata"),
                          type = c("free","series","grid"),
                          permute.strata = FALSE,
                          maxperm = 9999, minperm = 99,
                          mirror = FALSE, constant = FALSE,
                          ncol = NULL, nrow = NULL,
                          all.perms = NULL)
{
    if(missing(type))
        type <- "free"
    else
        type <- match.arg(type)
    out <- list(strata = strata, nperm = nperm, complete = complete,
                type = type, permute.strata = permute.strata,
                maxperm = maxperm, minperm = minperm,
                mirror = mirror, constant = constant,
                ncol = ncol, nrow = nrow, all.perms = all.perms,
                name.strata = deparse(substitute(strata)))
    class(out) <- "permControl"
    return(out)
}

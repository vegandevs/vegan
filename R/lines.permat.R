## S3 lines method for permat
`lines.permat` <-
function(x, type = "bray", ...)
{
    type <- match.arg(type, c("bray", "chisq"))
    if (type == "bray")
        toplot <- summary(x)$bray
    if (type == "chisq")
        toplot <- summary(x)$chisq$chisq.perm
    lines(toplot, ...)
}

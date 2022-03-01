## S3 lines method for permat
`lines.permat` <-
function(x, type = "bray", ...)
{
    lines(summary(x)[[match.arg(type, c("bray", "chisq"))]], ...)
}

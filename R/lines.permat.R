## S3 lines method for permat
`lines.permat` <-
function(x, ...)
{
    lines(summary(x)$bray, ...)
}

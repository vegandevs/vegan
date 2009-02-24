`lag.plot.permat` <-
function(x, type = "bray", ...)
{
    type <- match.arg(type, c("bray", "chisq"))
    ## this duality is required, because ylab can not be specified for lag.plot
    switch(type,
        "bray" = {
            out <- Bray.Curtis.Dissimilarities <- summary(x)[[type]]
            lag.plot(Bray.Curtis.Dissimilarities, ...)
        },
        "chisq" = {
            out <- Chi.Squared.Values <- summary(x)[[type]]
            lag.plot(Chi.Squared.Values, ...)
        })
#    lag.plot(summary(x)[[match.arg(type, c("bray", "chisq"))]], ...)
    invisible(out)
}

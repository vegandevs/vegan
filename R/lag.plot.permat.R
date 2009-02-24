`lag.plot.permat` <-
function(x, type = "bray", sub=NULL, ...)
{
    lags <- 1
    Box.type <- "Ljung-Box"
    type <- match.arg(type, c("bray", "chisq"))
    ## this duality is required, because ylab can not be specified for lag.plot
    switch(type,
        "bray" = {
            out <- Bray.Curtis.Dissimilarities <- summary(x)[[type]]
            if (is.null(sub)) {
                Bt <- Box.test(out, lag=lags, type=Box.type)[c(1,2,3)]
                names(Bt)[c(1,2)] <- c("X-squared", "df")
                Bt[[1]] <- round(Bt[[1]], 3)
                Bt[[3]] <- if (Bt[[3]] < 0.0001)
                    "< 0.0001" else round(Bt[[3]], 3)
                sub <- paste(Box.type, " test: ", paste(names(Bt), "=", Bt, collapse=", "), collapse="")
            }
            lag.plot(Bray.Curtis.Dissimilarities, lags=lags, sub=sub, ...)
            invisible(out)
        },
        "chisq" = {
            out <- Chi.Squared.Values <- summary(x)[[type]]
            if (is.null(sub)) {
                Bt <- Box.test(out, lag=lags, type=Box.type)[c(1,2,3)]
                names(Bt)[c(1,2)] <- c("X-squared", "df")
                Bt[[1]] <- round(Bt[[1]], 3)
                Bt[[3]] <- if (Bt[[3]] < 0.0001)
                    "< 0.0001" else round(Bt[[3]], 3)
                sub <- paste(Box.type, " test: ", paste(names(Bt), "=", Bt, collapse=", "), collapse="")
            }
            lag.plot(Chi.Squared.Values, lags=lags, sub=sub, ...)
            invisible(out)
        })
}

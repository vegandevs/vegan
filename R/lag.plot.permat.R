`lag.plot.permat` <-
function(x, type = "bray", sub=NULL, ...)
{
    type <- match.arg(type, c("bray", "chisq"))
    ## this duality is required, because ylab can not be specified for lag.plot
    switch(type,
        "bray" = {
            out <- Bray.Curtis.Dissimilarities <- summary(x)[[type]]
            if (is.null(sub)) {
                Bt <- Box.test(out)[c(1,2,3)]
                names(Bt)[c(1,2)] <- c("X-squared", "df")
                Bt[[1]] <- round(Bt[[1]], 3)
                Bt[[3]] <- if (Bt[[3]] < 0.0001)
                    "< 0.0001" else round(Bt[[3]], 3)
                sub <- paste("Box-Pierce test: ", paste(names(Bt), Bt, collapse=", "), collapse="")
            }
            lag.plot(Bray.Curtis.Dissimilarities, sub=sub, ...)
            invisible(out)
        },
        "chisq" = {
            out <- Chi.Squared.Values <- summary(x)[[type]]
            if (is.null(sub)) {
                Bt <- Box.test(out)[c(1,2,3)]
                names(Bt)[c(1,2)] <- c("X-squared", "df")
                Bt[[1]] <- round(Bt[[1]], 3)
                Bt[[3]] <- if (Bt[[3]] < 0.0001)
                    "< 0.0001" else round(Bt[[3]], 3)
                sub <- paste("Box-Pierce test: ", paste(names(Bt), Bt, collapse=", "), collapse="")
            }
            lag.plot(Chi.Squared.Values, sub=sub, ...)
            invisible(out)
        })
}

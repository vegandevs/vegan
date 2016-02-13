### Similar plotting functions as plot.prestonfit/fisherfit, but only
### for the data without the fitted models. These can be used with the
### result of as.preston(), as.fisher().

## as plot.prestonfit, but plots only the bars of as.preston()

`plot.preston` <-
        function (x, xlab = "Frequency", ylab = "Species", bar.col = "skyblue", 
                  ...) 
{
    freq <- x
    oct <- as.numeric(names(freq))
    noct <- max(oct) + 1
    plot(oct, freq, type = "n", ylim = c(0, max(freq)),
         xlim = c(-1, max(oct)), ylab = ylab, xlab = xlab, axes = FALSE, ...)
    axis(2)
    axis(1, at = 0:noct, labels = 2^(0:noct))
    box()
    rect(oct - 1, 0, oct, freq, col = bar.col, ...)
    invisible()
}


`lines.preston` <-
    function(x, xadjust = 0.5, ...)
{
    oct <- as.numeric(names(x)) - xadjust 
    lines(oct, x, ...)
}

## as plot.fisherfit, but plots only the bars of as.fisherfit

`plot.fisher` <-
    function(x, xlab = "Frequency", ylab = "Species", bar.col = "skyblue",
             kind = c("bar", "hiplot", "points", "lines"), add = FALSE,
             ...)
{
    kind <- match.arg(kind)
    freq <- as.numeric(names(x))
    if (!add)
        plot(freq, x, ylab=ylab, xlab=xlab,
             ylim=c(0,max(x)),  xlim=c(0.5, max(freq)+0.5), type="n", ...)
    switch(kind,
           "bar" = rect(freq-0.5,0,freq+0.5,x, col=bar.col, ...),
           "hiplot" = points(freq, x, col =bar.col, type = "h",  ...),
           "points" = points(freq, x, col =bar.col, ...),
           "lines" = lines(freq, x, col =bar.col, ...)
       )
    invisible()
}

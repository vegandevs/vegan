`plot.renyi` <-
    function(x, ...)
{
    ## one site: make into 1-row matrix
    if (!is.data.frame(x))
        x <- t(x)
    divs <- colnames(x)
    ndivs <- length(divs)
    matplot(seq_len(ndivs), t(x), axes=FALSE, type="l",
            col=palette(), lty=rep(1:6, each=length(palette())),
            xlab = "alpha", ylab = "Diversity", ...)
    axis(1, at = seq_len(ndivs), labels=divs)
    axis(2)
    box()
    if(nrow(x) > 1)
        legend("topright", rownames(x), lty=rep(1:6, each=length(palette())),
                                                col = palette())
}

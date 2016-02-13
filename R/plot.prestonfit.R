"plot.prestonfit" <-
    function (x, xlab = "Frequency", ylab = "Species", bar.col = "skyblue", 
              line.col = "red", lwd = 2, ...) 
{
    freq <- x$freq
    oct <- as.numeric(names(freq))
    noct <- max(oct) + 1
    plot(oct, freq, type = "n", ylim = c(0, max(freq)), xlim = c(-1, 
                                                        max(oct)), ylab = ylab, xlab = xlab, axes = FALSE, ...)
    axis(2)
    axis(1, at = 0:noct, labels = 2^(0:noct))
    box()
    rect(oct - 1, 0, oct, freq, col = bar.col, ...)
    p <- x$coefficients
    curve(p[3] * exp(-(x - p[1])^2/2/p[2]^2), -1, max(oct), add = TRUE, 
          col = line.col, lwd = lwd, ...)
    segments(p["mode"], 0, p["mode"], p["S0"], col = line.col, ...)
    segments(p["mode"] - p["width"], p["S0"] * exp(-0.5), p["mode"] + 
             p["width"], p["S0"] * exp(-0.5), col = line.col, ...)
    invisible()
}

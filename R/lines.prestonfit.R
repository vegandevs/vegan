"lines.prestonfit" <-
    function(x, line.col = "red", lwd = 2, ...)
{
    p <- x$coefficients
    freq <- x$freq
    oct <- as.numeric(names(freq))
    curve(p[3] * exp(-(x-p[1])^2/2/p[2]^2), -1, max(oct), add = TRUE,
          col = line.col, lwd = lwd, ...)
        segments(p["mode"], 0, p["mode"], p["S0"], col = line.col, 
        ...)
    segments(p["mode"] - p["width"], p["S0"] * exp(-0.5), p["mode"] + 
        p["width"], p["S0"] * exp(-0.5), col = line.col, ...)
    invisible()
}

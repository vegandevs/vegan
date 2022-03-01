`plot.betadiver` <-
    function (x, ...) 
{
    xy <- scores(x, ...)
    plot(c(0, 1), c(0, sqrt(0.75)), type = "n", axes = FALSE, 
         xlab = "", ylab = "", asp = 1)
    for (tic in seq(0.2, 0.8, by = 0.2)) {
        segments(tic, 0, tic/2, sqrt(0.75) * tic, lty = 3)
        segments(tic/2, sqrt(0.75) * tic, 1 - tic/2, sqrt(0.75) * 
                 tic, lty = 3)
        segments(tic, 0, tic/2 + 0.5, sqrt(0.75) * (1 - tic), 
                 lty = 3)
    }
    text(c(0, 1, 0.5), c(0, 0, sqrt(0.75)), c("b'", "c'", "a'"), 
         pos = c(2, 4, 3), cex = par("cex.axis"), xpd=TRUE)
    lines(c(0, 1, 0.5, 0), c(0, 0, sqrt(0.75), 0), xpd = TRUE)
    points(xy, ...)
    class(xy) <- "ordiplot"
    invisible(xy)
}

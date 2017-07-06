`coverscale` <-
    function (x, scale = c("Braun.Blanquet", "Domin", "Hult", "Hill",
                 "fix", "log"), maxabund, character = TRUE)
{
    scale <- match.arg(scale)
    sol <- as.data.frame(x)
    x <- as.matrix(x)
    switch(scale, Braun.Blanquet = {
        codes <- c("r", "+", as.character(1:5))
        lims <- c(0, 0.1, 1, 5, 25, 50, 75, 100)
    }, Domin = {
        codes <- c("+", as.character(1:9), "X")
        lims <- c(0, 0.01, 0.1, 1, 5, 10, 25, 33, 50, 75, 90,
                  100)
    }, Hult = {
        codes <- as.character(1:5)
        lims <- c(0, 100/2^(4:1), 100)
    }, Hill = {
        codes <- as.character(1:5)
        lims <- c(0, 2, 5, 10, 20, 100)
    }, fix = {
        codes <- c("+", as.character(1:9), "X")
        lims <- c(0:10, 11 - 10 * .Machine$double.eps)
    }, log = {
        codes <- c("+", as.character(1:9))
        if (missing(maxabund))
            maxabund <- max(x)
        lims <- c(0, maxabund/2^(9:1), maxabund)
    })
    for (i in 1:nrow(x)) {
        if (!character)
            codes <- FALSE
        tmp <- x[i, ] > 0
        sol[i, tmp] <- cut(x[i, tmp], breaks = lims, labels = codes,
                           right = FALSE, include.lowest = TRUE)
    }
    attr(sol, "scale") <-
        if (scale == "log")  paste("log, with maxabund", maxabund) else scale
    sol
}


"plot.radfit.frame" <-
    function (x, order.by, BIC = FALSE, model, legend = TRUE, as.table = TRUE, 
              ...) 
{
    modnam <- names(x[[1]]$models)
    if (!missing(model)) 
        pick <- pmatch(model, modnam, nomatch = FALSE)
    else pick <- FALSE
    pickmod <- function(x, pick, BIC) {
        if (pick) 
            return(pick)
        else {
            k <- if (BIC) 
                log(length(x$y))
            else 2
            which.min(AIC(x, k))
        }
    }
    Nhm <- length(x)
    Abundance <- unlist(lapply(x, function(x) x$y))
    Rank <- unlist(lapply(x, function(x) if (length(x$y) > 0) seq_along(x$y) else NULL))
    Site <- unlist(lapply(x, function(x) length(x$y)))
    N <- Site
    sitenames <- names(Site)
    Site <- rep(names(Site), Site)
    if (missing(order.by)) 
        order.by <- 1:Nhm
    else order.by <- order(order.by)
    Site <- factor(Site, levels = sitenames[order.by])
    fit <- unlist(lapply(x, function(x)
                         as.matrix(fitted(x))[, pickmod(x, 
                                                        pick, BIC)]))
    take <- sapply(x, function(x) pickmod(x, pick, BIC))
    take <- rep(take, N)
    cols <- trellis.par.get("superpose.line")$col
    cols <- cols[seq_along(cols)]
    if (legend) {
        mykey <- list(text = list(text = modnam), lines = list(lty = 1, 
                                                  col = cols[seq_along(modnam)], lwd = 2), columns = 3)
    }
    else {
        mykey <- NULL
    }
    tics <- function(x = max(Abundance), z = min(Abundance)) {
        ii <- round(c(log10(z), log10(x)))
        x10 <- 10^(ii[1]:ii[2])
        if (length(x10) < 3) 
            x10 <- c(outer(c(1, 2, 5), x10))
        else if (length(x10) < 6) 
            x10 <- c(outer(c(1, 3), x10))
        x10[x10 <= x & x10 >= z]
    }
    out <- xyplot(Abundance ~ Rank | Site, subscripts = TRUE, 
                  as.table = as.table, key = mykey, scales = list(y = list(log = 10, 
                                                                  at = tics())), panel = function(x, y, subscripts) {
                                                                      panel.xyplot(x, y, ...)
                                                                      panel.xyplot(x, log10(fit[subscripts]), type = "l", 
                                                                                   col = cols[take[min(subscripts)]], lwd = 2, ...)
                                                                  }, ...)
    out
}

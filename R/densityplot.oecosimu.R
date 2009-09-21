`densityplot.oecosimu` <-
    function(x, data, xlab = "Simulated", ...)
{
    require(lattice) || stop("requires package 'lattice'")
    sim <- t(x$oecosimu$simulated)
    obs <- x$oecosimu$statistic
    nm <- names(obs)[col(sim)]
    densityplot( ~ as.vector(sim) | factor(nm, levels = nm), xlab = xlab,
                panel = function(x, ...) {
                    panel.densityplot(x, ...)
                    panel.abline(v = obs[panel.number()], ...)
                },
                ...)
}

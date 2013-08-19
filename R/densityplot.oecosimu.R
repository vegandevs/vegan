`densityplot.oecosimu` <-
    function(x, data, xlab = "Simulated", ...)
{
    obs <- x$oecosimu$statistic
    sim <- rbind(obs, t(x$oecosimu$simulated))
    nm <- names(obs)[col(sim)]
    densityplot( ~ as.vector(sim) | factor(nm, levels = unique(nm)),
                xlab = xlab,
                panel = function(x, ...) {
                    panel.densityplot(x, ...)
                    panel.abline(v = obs[panel.number()], ...)
                },
                ...)
}
